import type { Feature, FeatureCollection } from 'geojson'
import MapLibreGL from 'maplibre-gl'
import 'maplibre-gl/dist/maplibre-gl.css'
import { useEffect, useRef, useState } from 'react'

interface MapViewProps {
  latitude: number
  longitude: number
  zoom?: number
  boundaryPolygon?: Feature | FeatureCollection | null
  rectangle?: { north: number; south: number; east: number; west: number } | null
  sidebarWidth?: number
  onMapClick?: (lat: number, lon: number) => void
  isSettingCenter?: boolean
  centerMarker?: { lat: number; lon: number } | null
  onMarkerDrag?: (lat: number, lon: number) => void
  boundaryBbox?: [number, number, number, number] | null
}

const DEFAULT_ZOOM = 13
const DEFAULT_STYLE = 'https://api.maptiler.com/maps/streets/style.json?key=HXTn4Pzqts6X4ceG00qZ'

const MapView: React.FC<MapViewProps> = ({
  latitude,
  longitude,
  zoom = DEFAULT_ZOOM,
  boundaryPolygon,
  rectangle,
  sidebarWidth,
  onMapClick,
  isSettingCenter,
  centerMarker,
  onMarkerDrag,
  boundaryBbox,
}) => {
  const mapRef = useRef<HTMLDivElement>(null)
  const mapInstance = useRef<MapLibreGL.Map | null>(null)
  const [markerPos, setMarkerPos] = useState<{ x: number; y: number } | null>(null)
  const markerRef = useRef<MapLibreGL.Marker | null>(null)
  const [isStyleLoaded, setIsStyleLoaded] = useState(false)

  useEffect(() => {
    if (!mapRef.current) return
    if (!mapInstance.current) {
      mapInstance.current = new MapLibreGL.Map({
        container: mapRef.current,
        style: DEFAULT_STYLE,
        center: [longitude, latitude],
        zoom,
      })
      mapInstance.current.on('load', () => {
        // No panBy here
      })
      mapInstance.current.on('style.load', () => {
        setIsStyleLoaded(true)
      })
    } else {
      mapInstance.current.setCenter([longitude, latitude])
      mapInstance.current.setZoom(zoom)
    }
    return () => {
      mapInstance.current?.remove()
      mapInstance.current = null
      setIsStyleLoaded(false)
    }
  }, [])

  useEffect(() => {
    if (mapInstance.current) {
      mapInstance.current.setCenter([longitude, latitude])
    }
  }, [latitude, longitude])

  useEffect(() => {
    if (mapInstance.current) {
      mapInstance.current.setZoom(zoom)
    }
  }, [zoom])

  // Add/Update boundary polygon layer
  useEffect(() => {
    if (!mapInstance.current || !boundaryPolygon || !isStyleLoaded) return
    const map = mapInstance.current
    if (map.getSource('boundary')) {
      ; (map.getSource('boundary') as unknown as MapLibreGL.GeoJSONSource).setData(boundaryPolygon)
    } else {
      map.addSource('boundary', {
        type: 'geojson',
        data: boundaryPolygon,
      })
      map.addLayer({
        id: 'boundary-fill',
        type: 'fill',
        source: 'boundary',
        paint: {
          'fill-color': '#d8dee9',
          'fill-opacity': 0.13,
        },
      })
      map.addLayer({
        id: 'boundary-line',
        type: 'line',
        source: 'boundary',
        paint: {
          'line-color': '#4c566a',
          'line-width': 3,
        },
      })
    }
    return () => {
      if (map && map.style && typeof map.getLayer === 'function' && map.getLayer('boundary-line'))
        map.removeLayer('boundary-line')
      if (map && map.style && typeof map.getLayer === 'function' && map.getLayer('boundary-fill'))
        map.removeLayer('boundary-fill')
      if (map && map.style && typeof map.getSource === 'function' && map.getSource('boundary'))
        map.removeSource('boundary')
    }
  }, [boundaryPolygon, isStyleLoaded])

  // Add/Update rectangle layer
  useEffect(() => {
    if (!mapInstance.current || !isStyleLoaded) return
    const map = mapInstance.current
    // Remove old rectangle layer/source if exists
    if (map.style && typeof map.getLayer === 'function' && map.getLayer('rectangle-layer'))
      map.removeLayer('rectangle-layer')
    if (map.style && typeof map.getSource === 'function' && map.getSource('rectangle'))
      map.removeSource('rectangle')
    if (!rectangle) return
    const { north, south, east, west } = rectangle
    const rectGeoJSON = {
      type: 'Feature' as const,
      geometry: {
        type: 'Polygon' as const,
        coordinates: [
          [
            [west, north],
            [east, north],
            [east, south],
            [west, south],
            [west, north],
          ],
        ],
      },
      properties: {},
    }
    map.addSource('rectangle', {
      type: 'geojson',
      data: rectGeoJSON,
    })
    map.addLayer({
      id: 'rectangle-layer',
      type: 'fill',
      source: 'rectangle',
      paint: {
        'fill-color': '#88c0d0',
        'fill-opacity': 0.18,
      },
    })
    map.addLayer({
      id: 'rectangle-outline',
      type: 'line',
      source: 'rectangle',
      paint: {
        'line-color': '#5e81ac',
        'line-width': 2.5,
      },
    })
    return () => {
      if (map.style && typeof map.getLayer === 'function' && map.getLayer('rectangle-layer'))
        map.removeLayer('rectangle-layer')
      if (map.style && typeof map.getLayer === 'function' && map.getLayer('rectangle-outline'))
        map.removeLayer('rectangle-outline')
      if (map.style && typeof map.getSource === 'function' && map.getSource('rectangle'))
        map.removeSource('rectangle')
    }
  }, [rectangle, isStyleLoaded])

  useEffect(() => {
    if (!mapInstance.current || !onMapClick || !isSettingCenter) return
    const map = mapInstance.current
    const handleClick = (e: unknown) => {
      const { lng, lat } = e.lngLat
      onMapClick(lat, lng)
    }
    map.on('click', handleClick)
    return () => {
      map.off('click', handleClick)
    }
  }, [onMapClick, isSettingCenter])

  useEffect(() => {
    if (!mapInstance.current || !centerMarker) {
      setMarkerPos(null)
      return
    }
    const map = mapInstance.current
    const { lat, lon } = centerMarker
    const { x, y } = map.project([lon, lat])
    setMarkerPos({ x, y })
  }, [centerMarker])

  useEffect(() => {
    if (!mapInstance.current) return
    const map = mapInstance.current
    // Remove old marker if exists
    if (markerRef.current) {
      markerRef.current.remove()
      markerRef.current = null
    }
    if (!centerMarker) return
    // Create a custom marker element
    const el = document.createElement('div')
    el.style.width = '32px'
    el.style.height = '32px'
    el.style.display = 'flex'
    el.style.alignItems = 'center'
    el.style.justifyContent = 'center'
    el.style.pointerEvents = 'auto'
    el.innerHTML = `<svg xmlns="http://www.w3.org/2000/svg" width="32" height="32" fill="none" viewBox="0 0 24 24"><path stroke="#5e81ac" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" d="M12 21c4.97-6.16 7.45-9.24 6.5-12.14C17.02 5.13 14.7 3 12 3s-5.02 2.13-6.5 5.86C4.55 11.76 7.03 14.84 12 21Z"/><circle cx="12" cy="10" r="3" stroke="#5e81ac" stroke-width="2"/></svg>`
    // Create marker
    const marker = new MapLibreGL.Marker({ element: el, anchor: 'center' })
      .setLngLat([centerMarker.lon, centerMarker.lat])
      .addTo(map)
    markerRef.current = marker
    return () => {
      marker.remove()
      markerRef.current = null
    }
  }, [centerMarker])

  useEffect(() => {
    if (!mapInstance.current || !boundaryBbox || !isStyleLoaded) return
    const map = mapInstance.current
    // [minLon, minLat, maxLon, maxLat]
    map.fitBounds(
      [
        [boundaryBbox[0], boundaryBbox[1]],
        [boundaryBbox[2], boundaryBbox[3]],
      ],
      {
        padding: 48,
        duration: 800,
        maxZoom: 15,
      },
    )
  }, [boundaryBbox, isStyleLoaded])

  useEffect(() => {
    // Promijeni cursor na map canvasu ovisno o isSettingCenter
    const canvas = document.querySelector('.maplibregl-canvas') as HTMLElement | null
    if (canvas) {
      if (isSettingCenter) {
        // Custom SVG crosshair, veći i jasno vidljiv
        canvas.style.cursor =
          "url(\"data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='48' height='48' viewBox='0 0 48 48'><circle cx='24' cy='24' r='8' fill='none' stroke='%235e81ac' stroke-width='2'/><line x1='24' y1='0' x2='24' y2='48' stroke='%235e81ac' stroke-width='2'/><line x1='0' y1='24' x2='48' y2='24' stroke='%235e81ac' stroke-width='2'/></svg>\") 24 24, crosshair"
      } else {
        canvas.style.cursor = 'grab'
      }
    }
    return () => {
      if (canvas) canvas.style.cursor = ''
    }
  }, [isSettingCenter])

  return (
    <div
      ref={mapRef}
      style={{ width: '100%', height: '100%', position: 'relative' }}
      id="map-view-container"
    />
  )
}

export default MapView
