"""
Patent model for storing patent information.
"""

from app import db


class Patent(db.Model):
    """Model for storing patent meta-information."""
    
    __tablename__ = 'patents'
    
    patent_id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    patent_number = db.Column(db.String(50), nullable=False, unique=True, index=True)
    title = db.Column(db.String(500), nullable=True)
    assignee = db.Column(db.String(300), nullable=True)
    filing_date = db.Column(db.Date, nullable=True)
    publication_date = db.Column(db.Date, nullable=True)
    abstract = db.Column(db.Text, nullable=True)
    created_at = db.Column(db.DateTime, default=db.func.current_timestamp())
    
    def to_dict(self):
        """Convert patent to dictionary representation."""
        return {
            'patent_id': self.patent_id,
            'patent_number': self.patent_number,
            'title': self.title,
            'assignee': self.assignee,
            'filing_date': self.filing_date.isoformat() if self.filing_date else None,
            'publication_date': self.publication_date.isoformat() if self.publication_date else None,
            'abstract': self.abstract,
            'created_at': self.created_at.isoformat() if self.created_at else None
        }
    
    def __repr__(self):
        return f'<Patent {self.patent_number}>'
