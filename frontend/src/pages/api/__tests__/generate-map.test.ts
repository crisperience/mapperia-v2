import { createMocks } from 'node-mocks-http'

import handler from '../generate-map'

describe('/api/generate-map', () => {
  it('should return 405 for GET', async () => {
    const { req, res } = createMocks({ method: 'GET' })
    await handler(req, res)
    expect(res.getStatusCode()).toBe(405)
  })

  it('should return 500 for backend error', async () => {
    const { req, res } = createMocks({
      method: 'POST',
      body: {
        lat: 0,
        lon: 0,
        width: 0,
        height: 0,
        formats: [],
        style: 'basic',
        location: '',
        layers: [],
      },
    })
    process.env.BACKEND_URL = 'http://localhost:9999' // force backend error
    await handler(req, res)
    expect(res.getStatusCode()).toBe(500)
  })

  // Note: For a real backend, you would mock fetch for a 200 test
})
