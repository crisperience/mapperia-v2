const nextConfig = {
  reactStrictMode: true,
  // Mapperia: route every mapbox-gl import to maplibre-gl
  webpack: config => {
    config.resolve.alias['mapbox-gl'] = 'maplibre-gl'
    return config
  },
  async rewrites() {
    return [
      {
        source: '/api/output/:path*',
        destination: 'http://localhost:8000/output/:path*', // Route to backend's static files
      },
      {
        source: '/api/:path*',
        destination: 'http://localhost:8000/api/:path*',
      },
    ]
  },
}

module.exports = nextConfig
