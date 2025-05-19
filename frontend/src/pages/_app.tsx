// Import maplibre-gl CSS for map styling
import 'maplibre-gl/dist/maplibre-gl.css'
import type { AppProps } from 'next/app'
import { Inter } from 'next/font/google'

// Simple global styles
import '../styles/globals.css'

const inter = Inter({ subsets: ['latin'], variable: '--font-inter' })

const App = ({ Component, pageProps }: AppProps) => (
  <main className={`${inter.variable} font-sans h-full`}>
    <Component {...pageProps} />
  </main>
)

export default App
