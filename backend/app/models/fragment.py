"""
Fragment model for storing RGroup replacement data.
"""

import json
from app import db


class Fragment(db.Model):
    """Model for storing molecular fragments."""
    
    __tablename__ = 'fragments'
    
    fragment_id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    smiles = db.Column(db.String(500), nullable=False, index=True)
    properties = db.Column(db.Text, nullable=True)  # JSON-encoded properties
    attachment_points = db.Column(db.Integer, default=1)
    name = db.Column(db.String(200), nullable=True)
    created_at = db.Column(db.DateTime, default=db.func.current_timestamp())
    updated_at = db.Column(db.DateTime, default=db.func.current_timestamp(), 
                          onupdate=db.func.current_timestamp())
    
    def get_properties(self):
        """Parse and return properties as a dictionary."""
        if self.properties:
            return json.loads(self.properties)
        return {}
    
    def set_properties(self, props_dict):
        """Set properties from a dictionary."""
        self.properties = json.dumps(props_dict)
    
    def to_dict(self):
        """Convert fragment to dictionary representation."""
        return {
            'fragment_id': self.fragment_id,
            'smiles': self.smiles,
            'properties': self.get_properties(),
            'attachment_points': self.attachment_points,
            'name': self.name,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'updated_at': self.updated_at.isoformat() if self.updated_at else None
        }
    
    def __repr__(self):
        return f'<Fragment {self.fragment_id}: {self.smiles}>'
