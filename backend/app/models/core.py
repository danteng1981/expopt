"""
Core model for storing scaffold/core structures.
"""

import json
from app import db


class Core(db.Model):
    """Model for storing molecular cores/scaffolds."""
    
    __tablename__ = 'cores'
    
    core_id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    smiles = db.Column(db.String(500), nullable=False, index=True)
    name = db.Column(db.String(200), nullable=True)
    description = db.Column(db.Text, nullable=True)
    properties = db.Column(db.Text, nullable=True)  # JSON-encoded properties
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
        """Convert core to dictionary representation."""
        return {
            'core_id': self.core_id,
            'smiles': self.smiles,
            'name': self.name,
            'description': self.description,
            'properties': self.get_properties(),
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'updated_at': self.updated_at.isoformat() if self.updated_at else None
        }
    
    def __repr__(self):
        return f'<Core {self.core_id}: {self.smiles}>'
