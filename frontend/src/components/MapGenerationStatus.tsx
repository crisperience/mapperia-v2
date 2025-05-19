import { Download } from 'lucide-react'
import React, { useEffect, useState } from 'react'

interface MapResult {
  svgUrl?: string
  pngUrl?: string
  layer_feature_counts?: Record<string, number>
  [key: string]: unknown
}

interface MapGenerationStatusProps {
  jobId: string
  onComplete?: (result: MapResult) => void
}

export const MapGenerationStatus: React.FC<MapGenerationStatusProps> = ({ jobId, onComplete }) => {
  const [status, setStatus] = useState<string>('processing')
  const [progress, setProgress] = useState<number>(0)
  const [error, setError] = useState<string | null>(null)
  const [result, setResult] = useState<MapResult | null>(null)

  useEffect(() => {
    const checkStatus = async () => {
      try {
        const response = await fetch(`/api/job-status/${jobId}`)
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`)
        }
        const data = await response.json()

        setStatus(data.status)
        setProgress(data.progress)

        if (data.error) {
          setError(data.error)
        }

        if (data.result) {
          setResult(data.result)
          if (onComplete) {
            onComplete(data.result)
          }
        }
      } catch (err) {
        if (process.env.NODE_ENV === 'development') {
          // eslint-disable-next-line no-console
          console.error('Error checking job status:', err)
        }
        setError(err instanceof Error ? err.message : 'Failed to check job status')
      }
    }

    const interval = setInterval(checkStatus, 1000)
    return () => clearInterval(interval)
  }, [jobId, onComplete])

  const handleDownload = (url: string) => {
    window.open(url, '_blank')
  }

  if (error) {
    return (
      <div className="mt-4 text-red-700 bg-red-100 border border-red-300 rounded p-4 text-center">
        {error}
      </div>
    )
  }

  return (
    <div className="mt-4 text-center">
      {status === 'processing' && (
        <div className="w-full flex flex-col items-center">
          <div className="relative w-16 h-16">
            <svg className="absolute top-0 left-0 w-full h-full" viewBox="0 0 64 64">
              <circle
                className="text-gray-200"
                strokeWidth="8"
                stroke="currentColor"
                fill="transparent"
                r="28"
                cx="32"
                cy="32"
              />
              <circle
                className="text-blue-500"
                strokeWidth="8"
                strokeDasharray={2 * Math.PI * 28}
                strokeDashoffset={2 * Math.PI * 28 * (1 - progress / 100)}
                strokeLinecap="round"
                stroke="currentColor"
                fill="transparent"
                r="28"
                cx="32"
                cy="32"
              />
            </svg>
            <span className="absolute inset-0 flex items-center justify-center text-lg font-semibold">
              {Math.round(progress)}%
            </span>
          </div>
          <div className="mt-2 text-gray-700">Generating map...</div>
        </div>
      )}

      {status === 'completed' && result && (
        <div>
          <div className="text-lg font-semibold mb-2">Map generated!</div>
          <div className="flex gap-4 justify-center mb-2">
            {result.svgUrl && (
              <button
                className="inline-flex items-center gap-2 px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 transition"
                onClick={() => handleDownload(result.svgUrl)}
                type="button"
                disabled={!result.svgUrl}
              >
                <Download size={18} /> Download SVG
              </button>
            )}
            {result.pngUrl && (
              <button
                className="inline-flex items-center gap-2 px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700 transition"
                onClick={() => handleDownload(result.pngUrl)}
                type="button"
                disabled={!result.pngUrl}
              >
                <Download size={18} /> Download PNG
              </button>
            )}
          </div>
          {result.layer_feature_counts && (
            <div className="mt-2">
              <div className="text-sm font-medium text-gray-700 mb-1">Layer Statistics:</div>
              {Object.entries(result.layer_feature_counts).map(([layer, count]) => (
                <div key={layer} className="text-xs text-gray-600">
                  {layer}: {count} features
                </div>
              ))}
            </div>
          )}
        </div>
      )}
    </div>
  )
}
