import type { Feature, FeatureCollection } from 'geojson'
import { Download, MapPin } from 'lucide-react'
import dynamic from 'next/dynamic'
import { Inter } from 'next/font/google'
import Head from 'next/head'
import Image from 'next/image'
import { useEffect, useRef, useState } from 'react'

import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { useSearchSuggestions } from '@/hooks/useSearchSuggestions'

const inter = Inter({ subsets: ['latin'], variable: '--font-inter' })

const MapView = dynamic(() => import('@/map/MapView'), { ssr: false })

// Types
type FormatKey = 'png' | 'svg_laser'
type Layer = { name: string; colorKey: string; color: string }
type LayersApiResponse = { layers: Layer[] }
type GenerateMapResponse = { job_id: string }
type JobStatusResponse = {
  status: 'completed' | 'failed' | string
  result?: { svgUrl?: string; pngUrl?: string }
  error?: string
  progress?: number
}

const LAYER_DESCRIPTIONS: Record<string, string> = {
  buildings: 'Buildings',
  greens: 'Greens',
  streets: 'Streets',
  blues: 'Blues',
  rail: 'Rail',
}

const LAYER_SUMMARIES: Record<string, string> = {
  buildings: 'houses, apartments',
  greens: 'parks, gardens',
  streets: 'roads, paths',
  blues: 'lakes, rivers',
  rail: 'railways, trams',
}

const LAYER_TOOLTIPS: Record<string, string> = {
  buildings: 'All OSM features with the building=* tag (houses, apartments, commercial, etc.)',
  greens:
    'Parks, gardens, golf courses, and recreation grounds (leisure=park|garden|golf_course|recreation_ground)',
  streets: 'All OSM features with the highway=* tag (roads, streets, paths, etc.)',
  blues: 'Water features (natural=water, waterway=*, landuse=reservoir, water=*)',
  rail: 'All OSM features with the railway=* tag (railways, trams, subways, etc.)',
}

const HomePage = () => {
  const inputRef = useRef<HTMLInputElement>(null)
  const [search, setSearch] = useState('')
  const [lat, setLat] = useState('')
  const [lon, setLon] = useState('')
  const [width, setWidth] = useState('0')
  const [height, setHeight] = useState('0')
  const [formats, setFormats] = useState<Record<FormatKey, boolean>>({
    png: false,
    svg_laser: false,
  })
  const [loading, setLoading] = useState(false)
  const [status, setStatus] = useState('')
  const [error, setError] = useState('')
  const [availableLayers, setAvailableLayers] = useState<Layer[]>([])
  const [layerConfig, setLayerConfig] = useState<Record<string, boolean>>({})
  const [downloadLinks, setDownloadLinks] = useState<{
    png?: string
    svgLaser?: string
  }>({})
  const [boundaryPolygon, setBoundaryPolygon] = useState<Feature | FeatureCollection | null>(null)
  const [latError, setLatError] = useState('')
  const [lonError, setLonError] = useState('')
  const [widthError, setWidthError] = useState('')
  const [heightError, setHeightError] = useState('')
  const [rectangle, setRectangle] = useState<{
    north: number
    south: number
    east: number
    west: number
  } | null>(null)
  const [isSettingCenter, setIsSettingCenter] = useState(false)
  const [centerMarker, setCenterMarker] = useState<{ lat: number; lon: number } | null>(null)
  const [boundaryBbox, setBoundaryBbox] = useState<[number, number, number, number] | null>(null)

  const {
    search: searchValue,
    setSearch: setSearchValue,
    suggestions: suggestionsHook,
    showSuggestions: showSuggestionsHook,
    searchLoading: searchLoadingHook,
    searchError: searchErrorHook,
    handleInputBlur: handleInputBlurHook,
    handleSuggestionClick: handleSuggestionClickHook,
    debouncedSearch: debouncedSearchHook,
  } = useSearchSuggestions({
    width,
    height,
    setLat,
    setLon,
    setBoundaryPolygon,
    setPreviewRectangle: undefined,
    setBoundaryBbox,
    inputRef,
  })

  useEffect(() => {
    fetch('/api/layers')
      .then(res => res.json())
      .then((data: LayersApiResponse) => {
        setAvailableLayers(data.layers)
        setLayerConfig(Object.fromEntries(data.layers.map(l => [l.name, false])))
      })
      .catch(() => setError('Failed to load layers'))
  }, [])

  const validateInputs = () => {
    let valid = true
    // Latitude
    const latNum = Number(lat)
    if (!lat || Number.isNaN(latNum) || latNum < -90 || latNum > 90) {
      setLatError('Latitude must be a number between -90 and 90')
      valid = false
    } else {
      setLatError('')
    }
    // Longitude
    const lonNum = Number(lon)
    if (!lon || Number.isNaN(lonNum) || lonNum < -180 || lonNum > 180) {
      setLonError('Longitude must be a number between -180 and 180')
      valid = false
    } else {
      setLonError('')
    }
    // Width (X)
    const widthNum = Number(width)
    if (!width || Number.isNaN(widthNum) || widthNum <= 0) {
      setWidthError('X (width) must be a positive number')
      valid = false
    } else {
      setWidthError('')
    }
    // Height (Y)
    const heightNum = Number(height)
    if (!height || Number.isNaN(heightNum) || heightNum <= 0) {
      setHeightError('Y (height) must be a positive number')
      valid = false
    } else {
      setHeightError('')
    }
    return valid
  }

  const handleGenerate = async () => {
    if (loading) return
    validateInputs()
    setLoading(true)
    setStatus('Generating map...')
    setError('')
    setDownloadLinks({})
    try {
      const res = await fetch('/api/generate-map', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          lat: parseFloat(lat),
          lon: parseFloat(lon),
          width: parseFloat(width),
          height: parseFloat(height),
          formats: Object.entries(formats)
            .filter(([_, v]) => v)
            .map(([k]) => k),
          style: 'basic',
          location: searchValue,
          layers: Object.keys(layerConfig).filter(k => layerConfig[k]),
        }),
      })
      if (!res.ok) throw new Error('Backend error')
      const data = await res.json()
      const newDownloadLinks = {
        png: data.pngUrl ? data.pngUrl : undefined,
        svgLaser: data.svgLaserUrl ? data.svgLaserUrl : undefined,
      }
      setDownloadLinks(newDownloadLinks)
      setStatus('Map generated!')
      setLoading(false)
    } catch (err) {
      setError('Failed to generate map')
      setStatus('')
      setLoading(false)
    }
  }

  // Sync rectangle with X/Y/center
  useEffect(() => {
    const latNum = Number(lat)
    const lonNum = Number(lon)
    const widthKm = Number(width)
    const heightKm = Number(height)
    if (
      !lat ||
      !lon ||
      !width ||
      !height ||
      Number.isNaN(latNum) ||
      Number.isNaN(lonNum) ||
      Number.isNaN(widthKm) ||
      Number.isNaN(heightKm) ||
      widthKm <= 0 ||
      heightKm <= 0
    ) {
      setRectangle(null)
      return
    }
    const dLat = heightKm / 2 / 111.32
    const dLon = widthKm / 2 / (111.32 * Math.cos((latNum * Math.PI) / 180))
    setRectangle({
      north: latNum + dLat,
      south: latNum - dLat,
      east: lonNum + dLon,
      west: lonNum - dLon,
    })
  }, [lat, lon, width, height])

  const handleMapClick = (clickedLat: number, clickedLon: number) => {
    setLat(clickedLat.toFixed(6))
    setLon(clickedLon.toFixed(6))
    setCenterMarker({ lat: clickedLat, lon: clickedLon })
  }

  const handleMarkerDrag = (newLat: number, newLon: number) => {
    setLat(newLat.toFixed(6))
    setLon(newLon.toFixed(6))
    setCenterMarker({ lat: newLat, lon: newLon })
  }

  return (
    <div
      className={`${inter.variable} font-sans text-[18px] relative w-screen h-screen overflow-hidden flex`}
    >
      <Head>
        <title>Mapperia - Map Generator</title>
      </Head>
      {/* Sidebar */}
      <aside className="h-full w-96 max-w-full bg-[var(--nord1)] text-[var(--nord6)] shadow-2xl border-r border-[var(--nord2)] z-10 flex flex-col">
        <div className="flex flex-col gap-10 p-12 h-full overflow-y-auto scrollbar-thin scrollbar-thumb-[var(--nord3)] scrollbar-track-[var(--nord0)]">
          {/* Logo */}
          <div className="flex justify-center items-center mb-4">
            <Image
              src="/logo.svg"
              alt="Mapperia Logo"
              width={240}
              height={72}
              className="mx-auto"
              priority
            />
          </div>
          <div className="flex-1 flex flex-col gap-8">
            {/* Location / Search */}
            <div className="flex flex-col relative">
              <Label htmlFor="search" className="text-[18px] font-medium text-[var(--nord6)] mb-1">
                Location
              </Label>
              <Input
                id="search"
                ref={inputRef}
                placeholder="Insert city or district..."
                value={searchValue}
                onChange={e => {
                  setSearchValue(e.target.value)
                  debouncedSearchHook(e.target.value)
                }}
                onBlur={handleInputBlurHook}
                className="mt-1 bg-[var(--nord2)] text-[var(--nord6)] border-[var(--nord3)] placeholder:text-[#eceff4] placeholder:opacity-70 placeholder:text-[16px] focus:border-[var(--nord6)] focus:ring-[var(--nord6)] rounded-[0.25rem]"
                autoComplete="off"
              />
              {lat &&
                lon &&
                !Number.isNaN(Number(lat)) &&
                !Number.isNaN(Number(lon)) &&
                Math.abs(Number(lat)) <= 90 &&
                Math.abs(Number(lon)) <= 180 && (
                  <div className="flex items-center gap-2 mt-4 select-none">
                    <span
                      role="button"
                      tabIndex={0}
                      className={`inline-flex items-center justify-center w-9 h-9 rounded-full transition-colors border-2 ${isSettingCenter
                        ? 'bg-[var(--nord10)] text-[var(--nord6)] border-[var(--nord9)] shadow-lg'
                        : 'bg-[var(--nord2)] text-[var(--nord6)] border-[var(--nord3)]'
                        } cursor-pointer`}
                      onClick={() => setIsSettingCenter(!isSettingCenter)}
                      onKeyDown={e =>
                        (e.key === 'Enter' || e.key === ' ') && setIsSettingCenter(!isSettingCenter)
                      }
                      aria-label="Set center on map"
                    >
                      <MapPin className="w-6 h-6" />
                    </span>
                    <span className="text-[16px] text-[var(--nord4)]">Set center on map</span>
                  </div>
                )}
              {searchLoadingHook && (
                <span className="text-xs text-[var(--nord4)] mt-1">Searching...</span>
              )}
              {searchErrorHook && (
                <span className="text-xs text-[var(--nord11)] mt-1">{searchErrorHook}</span>
              )}
              {showSuggestionsHook && (
                <div className="absolute top-full left-0 right-0 z-20 bg-[var(--nord2)] border border-[var(--nord3)] rounded-md shadow-lg max-h-56 overflow-auto mt-1">
                  {suggestionsHook.length > 0 ? (
                    suggestionsHook.map(s => (
                      <button
                        key={s.label}
                        type="button"
                        className="w-full text-left px-4 py-2 hover:bg-[var(--nord3)] text-[var(--nord6)] text-sm cursor-pointer"
                        onMouseDown={() => handleSuggestionClickHook(s)}
                      >
                        {s.label}
                      </button>
                    ))
                  ) : (
                    <div className="px-4 py-2 text-[var(--nord4)] text-sm">No results found</div>
                  )}
                </div>
              )}
            </div>
            <div className="border-t border-[var(--nord2)] my-2" />
            {/* Lat/Lon */}
            <div className="flex flex-col gap-4">
              <div>
                <Label htmlFor="lat" className="text-[18px] font-medium text-[var(--nord6)] mb-1">
                  Latitude (Y)
                </Label>
                <Input
                  id="lat"
                  type="text"
                  inputMode="decimal"
                  pattern="[0-9.\-]*"
                  placeholder="Insert latitude (e.g. 45.8150)"
                  value={lat}
                  onChange={e => setLat(e.target.value)}
                  className="mt-1 bg-[var(--nord2)] text-[var(--nord6)] text-[16px] border-[var(--nord3)] placeholder:text-[#eceff4] placeholder:opacity-70 placeholder:text-[16px] focus:border-[var(--nord6)] focus:ring-[var(--nord6)] rounded-[0.25rem]"
                  autoComplete="off"
                />
                {latError && <span className="text-xs text-[var(--nord11)] mt-1">{latError}</span>}
              </div>
              <div>
                <Label htmlFor="lon" className="text-[18px] font-medium text-[var(--nord6)] mb-1">
                  Longitude (X)
                </Label>
                <Input
                  id="lon"
                  type="text"
                  inputMode="decimal"
                  pattern="[0-9.\-]*"
                  placeholder="Insert longitude (e.g. 15.9819)"
                  value={lon}
                  onChange={e => setLon(e.target.value)}
                  className="mt-1 bg-[var(--nord2)] text-[var(--nord6)] text-[16px] border-[var(--nord3)] placeholder:text-[#eceff4] placeholder:opacity-70 placeholder:text-[16px] focus:border-[var(--nord6)] focus:ring-[var(--nord6)] rounded-[0.25rem]"
                  autoComplete="off"
                />
                {lonError && <span className="text-xs text-[var(--nord11)] mt-1">{lonError}</span>}
              </div>
            </div>
            <div className="border-t border-[var(--nord2)] my-2" />
            {/* X/Y Dimensions */}
            <div className="grid grid-cols-2 gap-6">
              <div>
                <Label htmlFor="width" className="text-[18px] font-medium text-[var(--nord6)] mb-1">
                  X (km)
                </Label>
                <Input
                  id="width"
                  type="number"
                  placeholder="0"
                  value={width}
                  onChange={e => setWidth(e.target.value)}
                  min={1}
                  step={1}
                  className="mt-1 bg-[var(--nord2)] text-[var(--nord6)] border-[var(--nord3)] placeholder:text-[#eceff4] placeholder:opacity-70 placeholder:text-[16px] focus:border-[var(--nord6)] focus:ring-[var(--nord6)] rounded-[0.25rem]"
                />
                {widthError && (
                  <span className="text-xs text-[var(--nord11)] mt-1">{widthError}</span>
                )}
              </div>
              <div>
                <Label
                  htmlFor="height"
                  className="text-[18px] font-medium text-[var(--nord6)] mb-1"
                >
                  Y (km)
                </Label>
                <Input
                  id="height"
                  type="number"
                  placeholder="0"
                  value={height}
                  onChange={e => setHeight(e.target.value)}
                  min={1}
                  step={1}
                  className="mt-1 bg-[var(--nord2)] text-[var(--nord6)] border-[var(--nord3)] placeholder:text-[#eceff4] placeholder:opacity-70 placeholder:text-[16px] focus:border-[var(--nord6)] focus:ring-[var(--nord6)] rounded-[0.25rem]"
                />
                {heightError && (
                  <span className="text-xs text-[var(--nord11)] mt-1">{heightError}</span>
                )}
              </div>
            </div>
            <div className="border-t border-[var(--nord2)] my-2" />
            {/* Layer Configuration */}
            <div className="flex flex-col gap-2">
              <span className="block text-[18px] font-medium text-[var(--nord6)] mb-1">
                Output Layers
              </span>
              <div className="flex flex-col gap-3">
                {availableLayers
                  .filter(layer => layer.name !== 'front_cover' && layer.name !== 'back_cover')
                  .map(layer => (
                    <label
                      key={layer.name}
                      htmlFor={`layer-${layer.name}`}
                      className="inline-flex items-center gap-2 cursor-pointer select-none"
                    >
                      <input
                        id={`layer-${layer.name}`}
                        type="checkbox"
                        checked={!!layerConfig[layer.name]}
                        onChange={e =>
                          setLayerConfig(prev => ({
                            ...prev,
                            [layer.name]: e.target.checked,
                          }))
                        }
                      />
                      <span
                        style={{
                          display: 'inline-block',
                          width: 18,
                          height: 18,
                          background: layer.color,
                          border: '1px solid #888',
                          marginRight: 8,
                          borderRadius: 3,
                        }}
                      />
                      <span className="text-[16px] text-[var(--nord6)]">
                        {LAYER_DESCRIPTIONS[layer.name] || layer.name}
                        <span className="text-[var(--nord4)] text-sm ml-2">
                          ({(LAYER_SUMMARIES[layer.name] || '').toLowerCase()})
                        </span>
                      </span>
                    </label>
                  ))}
              </div>
            </div>
            <div className="border-t border-[var(--nord2)] my-2" />
            {/* Output format */}
            <div className="flex flex-col gap-2">
              <span className="block text-[18px] font-medium text-[var(--nord6)] mb-2">
                Output Format
              </span>
              <div className="flex flex-col gap-2 w-full max-w-full">
                <label
                  htmlFor="format-png"
                  className="inline-flex items-center gap-2 cursor-pointer select-none"
                >
                  <input
                    id="format-png"
                    type="checkbox"
                    checked={formats.png}
                    onChange={e => setFormats(f => ({ ...f, png: e.target.checked }))}
                  />
                  <span className="text-[16px] text-[var(--nord6)]">
                    PNG{' '}
                    <span className="text-[var(--nord4)] text-sm ml-1">
                      (preview image, rasterized)
                    </span>
                  </span>
                </label>
                <label
                  htmlFor="format-svg-laser"
                  className="inline-flex items-center gap-2 cursor-pointer select-none"
                >
                  <input
                    id="format-svg-laser"
                    type="checkbox"
                    checked={formats.svg_laser}
                    onChange={e => setFormats(f => ({ ...f, svg_laser: e.target.checked }))}
                  />
                  <span className="text-[16px] text-[var(--nord6)]">
                    SVG{' '}
                    <span className="text-[var(--nord4)] text-sm ml-1">
                      (laser ready, vector, lines only)
                    </span>
                  </span>
                </label>
              </div>
            </div>
          </div>
          {/* Generate Button */}
          <div className="pt-4">
            <Button
              onClick={handleGenerate}
              disabled={loading}
              size="lg"
              className="w-full font-semibold rounded-[0.25rem] text-[18px] py-3"
            >
              {loading ? (
                <>
                  <svg
                    className="animate-spin -ml-1 mr-2 h-5 w-5 text-white"
                    xmlns="http://www.w3.org/2000/svg"
                    fill="none"
                    viewBox="0 0 24 24"
                  >
                    <circle
                      className="opacity-25"
                      cx="12"
                      cy="12"
                      r="10"
                      stroke="currentColor"
                      strokeWidth="4"
                    />
                    <path
                      className="opacity-75"
                      fill="currentColor"
                      d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                    />
                  </svg>
                  GENERATE MAP
                </>
              ) : (
                'GENERATE MAP'
              )}
            </Button>
            {status && (
              <div className="mt-2 text-base text-[var(--nord6)] bg-[var(--nord2)] p-2 rounded-md border border-[var(--nord3)] flex flex-col items-center justify-center text-center w-full">
                {status}
                <div className="flex flex-row flex-wrap gap-2 mt-2 justify-center items-center">
                  {downloadLinks.png && (
                    <Button
                      asChild
                      className="w-full bg-green-600 hover:bg-green-700 text-white rounded-[0.25rem] px-3 py-2 text-sm"
                      variant="default"
                    >
                      <a
                        href={
                          downloadLinks.png
                            ? `${process.env.NEXT_PUBLIC_BACKEND_URL || 'http://localhost:8000'
                            }/api/output/${encodeURIComponent(
                              downloadLinks.png.split('/').pop() ?? '',
                            )}`
                            : ''
                        }
                        download
                      >
                        <Download className="mr-2 h-4 w-4" />
                        Download PNG
                      </a>
                    </Button>
                  )}
                  {downloadLinks.svgLaser && (
                    <Button
                      asChild
                      className="w-full bg-purple-600 hover:bg-purple-700 text-white rounded-[0.25rem] px-3 py-2 text-sm"
                      variant="default"
                    >
                      <a
                        href={
                          downloadLinks.svgLaser
                            ? `${process.env.NEXT_PUBLIC_BACKEND_URL || 'http://localhost:8000'
                            }/api/output/${encodeURIComponent(
                              downloadLinks.svgLaser.split('/').pop() ?? '',
                            )}`
                            : ''
                        }
                        download
                      >
                        <Download className="mr-2 h-4 w-4" />
                        Download SVG
                      </a>
                    </Button>
                  )}
                </div>
              </div>
            )}
            {error && (
              <div className="mt-4 text-base text-[var(--nord11)] bg-[var(--nord1)] p-3 rounded-md border border-[var(--nord11)]">
                {error}
              </div>
            )}
          </div>
        </div>
      </aside>
      {/* Map background */}
      <div className="flex-1 h-full z-0 relative">
        <MapView
          latitude={
            Number.isFinite(Number(lat)) && Math.abs(Number(lat)) <= 90 ? Number(lat) : 45.815
          }
          longitude={
            Number.isFinite(Number(lon)) && Math.abs(Number(lon)) <= 180 ? Number(lon) : 15.9819
          }
          boundaryPolygon={boundaryPolygon}
          rectangle={rectangle}
          onMapClick={isSettingCenter ? handleMapClick : undefined}
          isSettingCenter={isSettingCenter}
          centerMarker={isSettingCenter && centerMarker ? centerMarker : null}
          onMarkerDrag={handleMarkerDrag}
          sidebarWidth={384}
          boundaryBbox={boundaryBbox}
        />
        {!lat && !lon && (
          <div className="absolute inset-0 z-20 flex items-center justify-center bg-[var(--nord0)] bg-opacity-80">
            <div className="text-center">
              <h1 className="text-3xl font-bold text-[var(--nord6)] mb-4">Welcome to Mapperia!</h1>
              <p className="text-lg text-[var(--nord4)]">
                Insert a location or latitude and longitude to begin.
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}

export default HomePage
