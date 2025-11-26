"""
Routes package for ExpOpt API.
"""

from app.routes.fragments import fragments_bp
from app.routes.cores import cores_bp
from app.routes.tasks import tasks_bp

__all__ = ['fragments_bp', 'cores_bp', 'tasks_bp']
