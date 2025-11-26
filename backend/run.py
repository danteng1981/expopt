"""
ExpOpt Backend Application Entry Point

Run this file to start the Flask development server.
"""

import os
from app import create_app

app = create_app()

if __name__ == '__main__':
    # Only enable debug mode when explicitly set in environment
    debug_mode = os.environ.get('FLASK_DEBUG', 'false').lower() == 'true'
    app.run(debug=debug_mode, host='0.0.0.0', port=5000)
