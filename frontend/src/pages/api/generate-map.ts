import type { NextApiRequest, NextApiResponse } from 'next'
import fetch from 'node-fetch'

// Use env variable or fallback to localhost:8000
const BACKEND_URL = process.env.BACKEND_URL || 'http://localhost:8000'

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  if (req.method !== 'POST') {
    res.setHeader('Allow', ['POST'])
    return res.status(405).json({ error: 'Method Not Allowed' })
  }

  try {
    // Transform the request to match backend's expected format
    const backendPayload = {
      lat: req.body.lat,
      lon: req.body.lon,
      width: req.body.width,
      height: req.body.height,
      formats: req.body.formats,
      style: req.body.style,
      location: req.body.location,
      layers: req.body.layers,
    }

    try {
      const backendRes = await fetch(`${BACKEND_URL}/api/generate-map`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(backendPayload),
        timeout: 0, // no timeout
      })

      if (!backendRes.ok) {
        return res.status(500).json({ error: `Backend responded with status ${backendRes.status}` })
      }
      const data = await backendRes.json()
      return res.status(200).json(data)
    } catch (error: unknown) {
      return res.status(500).json({
        error: 'Failed to connect to backend',
        detail: error instanceof Error ? error.message : error,
      })
    }
  } catch (error: unknown) {
    return res.status(500).json({
      error: 'Failed to process request',
      detail: error instanceof Error ? error.message : error,
    })
  }
}
