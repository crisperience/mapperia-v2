/* eslint-disable consistent-return */
import fs from 'fs'
import { NextApiRequest, NextApiResponse } from 'next'
import path from 'path'

export default function handler(req: NextApiRequest, res: NextApiResponse) {
  const { output: filePathArr } = req.query
  if (!filePathArr || !Array.isArray(filePathArr)) {
    return res.status(400).json({ error: 'Missing file path' })
  }
  const fileName = filePathArr.join('/')
  const backendOutputPath = path.join(process.cwd(), '..', 'backend', 'output', fileName)
  if (process.env.NODE_ENV === 'development') {
    // eslint-disable-next-line no-console
    console.log('[DEBUG] Looking for file at:', backendOutputPath)
  }

  if (!fs.existsSync(backendOutputPath)) {
    if (process.env.NODE_ENV === 'development') {
      // eslint-disable-next-line no-console
      console.log('[DEBUG] File not found at:', backendOutputPath)
    }
    return res.status(404).json({ error: 'File not found' })
  }

  const ext = path.extname(fileName).toLowerCase()
  let contentType = 'application/octet-stream'
  if (ext === '.svg') contentType = 'image/svg+xml'
  if (ext === '.png') contentType = 'image/png'

  res.setHeader('Content-Type', contentType)
  res.setHeader('Content-Disposition', `attachment; filename="${path.basename(fileName)}"`)
  res.setHeader('Cache-Control', 'no-cache, no-store, must-revalidate')
  res.setHeader('Pragma', 'no-cache')
  res.setHeader('Expires', '0')

  try {
    const fileStream = fs.createReadStream(backendOutputPath)
    fileStream.pipe(res)
  } catch (error) {
    if (process.env.NODE_ENV === 'development') {
      // eslint-disable-next-line no-console
      console.error('[DEBUG] Error reading file:', error)
    }
    return res.status(500).json({ error: 'Error reading file' })
  }
}
