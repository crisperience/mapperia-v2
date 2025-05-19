import type { NextApiRequest, NextApiResponse } from 'next'

const BACKEND_URL = process.env.BACKEND_URL || 'http://localhost:8000'

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  if (req.method !== 'GET') {
    res.setHeader('Allow', ['GET'])
    return res.status(405).json({ error: 'Method Not Allowed' })
  }

  const { jobId } = req.query

  if (!jobId || typeof jobId !== 'string') {
    return res.status(400).json({ error: 'Invalid job ID' })
  }

  try {
    const controller = new AbortController()
    const timeout = setTimeout(() => controller.abort(), 120000) // 2 minutes

    const backendRes = await fetch(`${BACKEND_URL}/api/job-status/${jobId}`, {
      signal: controller.signal,
    })
    clearTimeout(timeout)

    if (!backendRes.ok) {
      return res.status(500).json({ error: `Backend responded with status ${backendRes.status}` })
    }
    const data = await backendRes.json()
    return res.status(200).json(data)
  } catch (error: unknown) {
    if (process.env.NODE_ENV === 'development') {
      // eslint-disable-next-line no-console
      console.error('Error fetching job status:', error)
    }
    return res.status(500).json({
      error: 'Failed to fetch job status',
      detail: error instanceof Error ? error.message : error,
    })
  }
}
