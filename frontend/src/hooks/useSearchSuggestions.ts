import bbox from '@turf/bbox'
import type { Feature, FeatureCollection, Geometry } from 'geojson'
import debounce from 'lodash.debounce'
import { useMemo, useState } from 'react'

interface Suggestion {
  label: string
  lat: string
  lon: string
  osm_type?: string
  osm_id?: string
}

interface UseSearchSuggestionsProps {
  width: string
  height: string
  setLat: (lat: string) => void
  setLon: (lon: string) => void
  setBoundaryPolygon: (poly: Feature | FeatureCollection | null) => void
  setPreviewRectangle?: (
    rect: { north: number; south: number; east: number; west: number } | null,
  ) => void
  inputRef: React.RefObject<HTMLInputElement>
}

export function useSearchSuggestions({
  width,
  height,
  setLat,
  setLon,
  setBoundaryPolygon,
  setPreviewRectangle = () => {},
  setBoundaryBbox = () => {},
  inputRef,
}: UseSearchSuggestionsProps & {
  setBoundaryBbox?: (bbox: [number, number, number, number] | null) => void
}) {
  const [search, setSearch] = useState('')
  const [suggestions, setSuggestions] = useState<Suggestion[]>([])
  const [showSuggestions, setShowSuggestions] = useState(false)
  const [searchLoading, setSearchLoading] = useState(false)
  const [searchError, setSearchError] = useState('')

  const debouncedSearch = useMemo(
    () =>
      debounce(async (value: string) => {
        setSearchError('')
        if (value.length < 3) {
          setSuggestions([])
          setShowSuggestions(false)
          return
        }
        setSearchLoading(true)
        try {
          const results = await fetch(
            `https://nominatim.openstreetmap.org/search?q=${encodeURIComponent(
              value,
            )}&format=json&addressdetails=1&extratags=1&limit=5`,
          ).then(res => res.json())
          const filtered = results.filter((s: Record<string, unknown>) => {
            const label = typeof s.display_name === 'string' ? s.display_name.toLowerCase() : ''
            const query = value.trim().toLowerCase()
            return label.split(/[^a-zA-Z0-9šđčćž]+/).some(word => word === query)
          })
          const mapped = filtered.map((s: Record<string, unknown>) => ({
            label: typeof s.display_name === 'string' ? s.display_name : '',
            lat: s.lat as string,
            lon: s.lon as string,
            osm_type: s.osm_type as string,
            osm_id: s.osm_id as string,
          }))
          setSuggestions(mapped)
          setShowSuggestions(mapped.length > 0)
        } catch (err) {
          setSearchError('Error searching location.')
          setSuggestions([])
          setShowSuggestions(false)
        } finally {
          setSearchLoading(false)
        }
      }, 600),
    [],
  )

  const handleSuggestionClick = async (s: Suggestion) => {
    setSearch(s.label)
    setLat(s.lat)
    setLon(s.lon)
    setShowSuggestions(false)
    setSuggestions([])
    inputRef.current?.blur()
    let foundPolygon = false
    try {
      let overpassQuery = ''
      if (s.osm_type && s.osm_id) {
        if (s.osm_type === 'relation') {
          overpassQuery = `[out:json][timeout:25];relation(${s.osm_id});out body;>;out skel qt;`
        } else if (s.osm_type === 'way') {
          overpassQuery = `[out:json][timeout:25];way(${s.osm_id});out body;>;out skel qt;`
        } else if (s.osm_type === 'node') {
          overpassQuery = `[out:json][timeout:25];node(${s.osm_id});out body;>;out skel qt;`
        }
      } else {
        overpassQuery = `[out:json];relation["name"="${s.label}"]["boundary"="administrative"];out geom;`
      }
      const url = `https://overpass-api.de/api/interpreter?data=${encodeURIComponent(
        overpassQuery,
      )}`
      const res = await fetch(url)
      const data = await res.json()
      if (data.elements && data.elements.length > 0) {
        const coords = data.elements[0].members
          .filter((m: Record<string, unknown>) => Array.isArray(m.geometry))
          .map((m: Record<string, unknown>) =>
            (m.geometry as Array<{ lon: number; lat: number }>).map(g => [g.lon, g.lat]),
          )
        let polygon
        if (coords.length === 1) {
          polygon = {
            type: 'Feature' as const,
            geometry: {
              type: 'Polygon' as const,
              coordinates: coords,
            },
            properties: {},
          }
        } else {
          polygon = {
            type: 'Feature' as const,
            geometry: {
              type: 'MultiPolygon' as const,
              coordinates: coords.map(c => [c]),
            },
            properties: {},
          }
        }
        if (
          polygon &&
          (polygon.geometry.type === 'Polygon' || polygon.geometry.type === 'MultiPolygon') &&
          Array.isArray(polygon.geometry.coordinates) &&
          polygon.geometry.coordinates.length > 0
        ) {
          setBoundaryPolygon(polygon)
          const bounds = bbox(polygon)
          if (
            Array.isArray(bounds) &&
            bounds.length === 4 &&
            bounds.every(v => Number.isFinite(v))
          ) {
            setBoundaryBbox(bounds as [number, number, number, number])
          } else {
            setBoundaryBbox(null)
          }
          foundPolygon = true
        } else {
          setBoundaryPolygon(null)
          setBoundaryBbox(null)
        }
      }
    } catch {}
    // Fallback to Nominatim if Overpass fails
    if (!foundPolygon) {
      try {
        const nominatimUrl = `https://nominatim.openstreetmap.org/search?format=json&polygon_geojson=1&q=${encodeURIComponent(
          s.label,
        )}`
        const nominatimRes = await fetch(nominatimUrl)
        const nominatimData = await nominatimRes.json()
        if (nominatimData[0]?.geojson) {
          const { geojson } = nominatimData[0]
          const feature: Feature = {
            type: 'Feature',
            geometry: geojson as Geometry,
            properties: {},
          }
          if (feature.geometry.type === 'Polygon' || feature.geometry.type === 'MultiPolygon') {
            if (
              Array.isArray(feature.geometry.coordinates) &&
              feature.geometry.coordinates.length > 0
            ) {
              setBoundaryPolygon(feature)
              const bounds = bbox(feature)
              if (
                Array.isArray(bounds) &&
                bounds.length === 4 &&
                bounds.every(v => Number.isFinite(v))
              ) {
                setBoundaryBbox(bounds as [number, number, number, number])
              } else {
                setBoundaryBbox(null)
              }
            } else {
              setBoundaryPolygon(null)
              setBoundaryBbox(null)
            }
          } else {
            setBoundaryPolygon(null)
            setBoundaryBbox(null)
          }
        } else {
          setBoundaryPolygon(null)
          if (typeof setBoundaryBbox === 'function') setBoundaryBbox(null)
        }
      } catch {
        setBoundaryPolygon(null)
        if (typeof setBoundaryBbox === 'function') setBoundaryBbox(null)
      }
    }
    // Calculate preview rectangle
    const latNum = parseFloat(s.lat)
    const lonNum = parseFloat(s.lon)
    const widthKm = parseFloat(width)
    const heightKm = parseFloat(height)
    // 1 deg latitude ~ 110.574 km, 1 deg longitude ~ 111.320*cos(lat) km
    const dLat = heightKm / 2 / 110.574
    const dLon = widthKm / 2 / (111.32 * Math.cos((latNum * Math.PI) / 180))
    setPreviewRectangle({
      north: latNum + dLat,
      south: latNum - dLat,
      east: lonNum + dLon,
      west: lonNum - dLon,
    })
  }

  const handleInputBlur = () => {
    setTimeout(() => setShowSuggestions(false), 100)
  }

  return {
    search,
    setSearch,
    suggestions,
    showSuggestions,
    searchLoading,
    searchError,
    handleInputBlur,
    handleSuggestionClick,
    debouncedSearch,
  }
}
