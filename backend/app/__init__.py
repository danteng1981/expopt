"""
ExpOpt Backend Application

This module initializes the Flask application with database support
for RGroup replacement and CoreHopping data.
"""

from flask import Flask
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()


def create_app(config_name=None):
    """Create and configure the Flask application."""
    app = Flask(__name__)
    
    # Configure the application
    if config_name == 'testing':
        app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///:memory:'
        app.config['TESTING'] = True
    else:
        app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///expopt.db'
    
    app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
    app.config['SECRET_KEY'] = 'expopt-secret-key'
    
    # Initialize extensions
    CORS(app)
    db.init_app(app)
    
    # Register blueprints
    from app.routes.fragments import fragments_bp
    from app.routes.cores import cores_bp
    from app.routes.tasks import tasks_bp
    
    app.register_blueprint(fragments_bp, url_prefix='/api/fragments')
    app.register_blueprint(cores_bp, url_prefix='/api/cores')
    app.register_blueprint(tasks_bp, url_prefix='/api/tasks')
    
    # Create database tables
    with app.app_context():
        db.create_all()
    
    return app
