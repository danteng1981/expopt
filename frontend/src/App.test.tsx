import React from 'react';
import { render, screen } from '@testing-library/react';
import App from './App';

test('renders ExpOpt navigation brand', () => {
  render(<App />);
  const brandElements = screen.getAllByText(/ExpOpt/i);
  expect(brandElements.length).toBeGreaterThan(0);
});
