"""
Task models for storing replacement/hopping task information and results.
"""

import json
from app import db


class Task(db.Model):
    """Model for storing RGroup/CoreHopping task submissions."""
    
    __tablename__ = 'tasks'
    
    task_id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    task_type = db.Column(db.String(50), nullable=False)  # 'rgroup' or 'core_hopping'
    input_smiles = db.Column(db.Text, nullable=False)
    constraints = db.Column(db.Text, nullable=True)  # JSON-encoded constraints
    status = db.Column(db.String(20), default='pending')  # pending, running, completed, failed
    created_at = db.Column(db.DateTime, default=db.func.current_timestamp())
    completed_at = db.Column(db.DateTime, nullable=True)
    error_message = db.Column(db.Text, nullable=True)
    
    # Relationship to results
    results = db.relationship('TaskResult', backref='task', lazy=True, cascade='all, delete-orphan')
    
    def get_constraints(self):
        """Parse and return constraints as a dictionary."""
        if self.constraints:
            return json.loads(self.constraints)
        return {}
    
    def set_constraints(self, constraints_dict):
        """Set constraints from a dictionary."""
        self.constraints = json.dumps(constraints_dict)
    
    def to_dict(self, include_results=False):
        """Convert task to dictionary representation."""
        data = {
            'task_id': self.task_id,
            'task_type': self.task_type,
            'input_smiles': self.input_smiles,
            'constraints': self.get_constraints(),
            'status': self.status,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'completed_at': self.completed_at.isoformat() if self.completed_at else None,
            'error_message': self.error_message
        }
        if include_results:
            data['results'] = [r.to_dict() for r in self.results]
        return data
    
    def __repr__(self):
        return f'<Task {self.task_id}: {self.task_type}>'


class TaskResult(db.Model):
    """Model for storing task results."""
    
    __tablename__ = 'task_results'
    
    result_id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    task_id = db.Column(db.Integer, db.ForeignKey('tasks.task_id'), nullable=False)
    output_smiles = db.Column(db.String(500), nullable=False)
    score = db.Column(db.Float, nullable=True)
    properties = db.Column(db.Text, nullable=True)  # JSON-encoded properties
    rank = db.Column(db.Integer, nullable=True)
    
    def get_properties(self):
        """Parse and return properties as a dictionary."""
        if self.properties:
            return json.loads(self.properties)
        return {}
    
    def set_properties(self, props_dict):
        """Set properties from a dictionary."""
        self.properties = json.dumps(props_dict)
    
    def to_dict(self):
        """Convert result to dictionary representation."""
        return {
            'result_id': self.result_id,
            'task_id': self.task_id,
            'output_smiles': self.output_smiles,
            'score': self.score,
            'properties': self.get_properties(),
            'rank': self.rank
        }
    
    def __repr__(self):
        return f'<TaskResult {self.result_id} for Task {self.task_id}>'
