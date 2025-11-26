/**
 * Molecule viewer component for displaying molecular structures
 * Uses SVG placeholder - in production, integrate RDKit.js or similar
 */

import React from 'react';
import './MoleculeViewer.css';

interface MoleculeViewerProps {
  smiles: string;
  width?: number;
  height?: number;
}

const MoleculeViewer: React.FC<MoleculeViewerProps> = ({ 
  smiles, 
  width = 200, 
  height = 150 
}) => {
  // Simple placeholder visualization
  // In production, use RDKit.js or a molecular rendering library
  return (
    <div className="molecule-viewer" style={{ width, height }}>
      <div className="molecule-placeholder">
        <svg viewBox="0 0 100 100" width={width * 0.8} height={height * 0.7}>
          {/* Simple hexagon to represent a ring structure */}
          <polygon 
            points="50,10 85,30 85,70 50,90 15,70 15,30" 
            fill="none" 
            stroke="#3498db" 
            strokeWidth="2"
          />
          {/* Atoms represented as circles */}
          <circle cx="50" cy="10" r="5" fill="#e74c3c" />
          <circle cx="85" cy="30" r="5" fill="#2ecc71" />
          <circle cx="85" cy="70" r="5" fill="#f39c12" />
          <circle cx="50" cy="90" r="5" fill="#9b59b6" />
          <circle cx="15" cy="70" r="5" fill="#1abc9c" />
          <circle cx="15" cy="30" r="5" fill="#e67e22" />
        </svg>
        <div className="smiles-label" title={smiles}>
          {smiles.length > 20 ? `${smiles.substring(0, 20)}...` : smiles}
        </div>
      </div>
    </div>
  );
};

export default MoleculeViewer;
