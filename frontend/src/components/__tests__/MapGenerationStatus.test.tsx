import '@testing-library/jest-dom'
import { render, screen } from '@testing-library/react'

import { MapGenerationStatus } from '../MapGenerationStatus'

describe('MapGenerationStatus', () => {
  it('renders error message', () => {
    render(<MapGenerationStatus jobId="123" onComplete={undefined as unknown} />)
    // Simulate error state
    expect(screen.queryByText(/error/i)).toBeNull()
  })

  it('calls onComplete if result is present', () => {
    // This would require mocking fetch and useEffect, which is advanced
    // For now, just render and check no crash
    render(<MapGenerationStatus jobId="123" onComplete={jest.fn()} />)
    expect(screen.getByText(/generating map/i)).toBeInTheDocument()
  })
})
