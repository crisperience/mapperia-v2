import type { NextApiRequest, NextApiResponse } from 'next'

interface Layer {
  name: string
  lcolor: [string, string]
}

interface LayerResponse {
  name: string
  colorKey: string
  color: string
}

interface ApiResponse {
  layers: LayerResponse[]
  colors: Record<string, string>
}

// Official hardcoded list of layers and colors (from backend/card.json and backend/lightburncolor.json)
const layers: Layer[] = [
  { name: 'buildings', lcolor: ['02red', '#FF0000'] },
  { name: 'greens', lcolor: ['03green', '#00FF00'] },
  { name: 'streets', lcolor: ['08lightgrey', '#B4B4B4'] },
  { name: 'blues', lcolor: ['01blue', '#0000FF'] },
  { name: 'rail', lcolor: ['09darkgrey', '#808080'] },
  { name: 'front_cover', lcolor: ['20darkred', '#D33F6A'] },
  { name: 'back_cover', lcolor: ['20darkred', '#D33F6A'] },
]

const colors: Record<string, string> = {
  '00black': '#000000',
  '01blue': '#0000FF',
  '02red': '#FF0000',
  '03green': '#00E000',
  '05orange': '#FF8000',
  '08lightgrey': '#B4B4B4',
  '15purple': '#A000A0',
  '16grey': '#808080',
  '17darkblue': '#7D87B9',
  '18brown': '#BB7784',
  '20darkred': '#D33F6A',
  '09darkgrey': '#808080',
}

export default function handler(
  req: NextApiRequest,
  res: NextApiResponse<ApiResponse | { error: string }>,
) {
  res.status(200).json({
    layers: layers.map(layer => ({
      name: layer.name,
      colorKey: layer.lcolor[0],
      color: layer.lcolor[1],
    })),
    colors,
  })
}
