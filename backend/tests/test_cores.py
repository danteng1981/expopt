"""
Tests for Core API endpoints.
"""

import json
import pytest


class TestCoresAPI:
    """Test cases for cores API."""
    
    def test_get_cores_empty(self, client):
        """Test getting cores when database is empty."""
        response = client.get('/api/cores')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['cores'] == []
        assert data['total'] == 0
    
    def test_create_core(self, client):
        """Test creating a new core."""
        core_data = {
            'smiles': 'c1ccccc1',
            'name': 'Benzene Core',
            'description': 'Simple benzene scaffold',
            'properties': {'MW': 78.11}
        }
        response = client.post(
            '/api/cores',
            data=json.dumps(core_data),
            content_type='application/json'
        )
        assert response.status_code == 201
        data = json.loads(response.data)
        assert data['smiles'] == 'c1ccccc1'
        assert data['name'] == 'Benzene Core'
    
    def test_create_core_missing_smiles(self, client):
        """Test creating core without SMILES returns error."""
        response = client.post(
            '/api/cores',
            data=json.dumps({'name': 'Test'}),
            content_type='application/json'
        )
        assert response.status_code == 400
    
    def test_get_core_by_id(self, client):
        """Test getting a specific core."""
        # Create core first
        core_data = {'smiles': 'c1ccncc1', 'name': 'Pyridine'}
        create_response = client.post(
            '/api/cores',
            data=json.dumps(core_data),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Get by ID
        response = client.get(f'/api/cores/{created["core_id"]}')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['smiles'] == 'c1ccncc1'
    
    def test_delete_core(self, client):
        """Test deleting a core."""
        # Create core first
        create_response = client.post(
            '/api/cores',
            data=json.dumps({'smiles': 'c1ccc2ccccc2c1'}),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Delete
        response = client.delete(f'/api/cores/{created["core_id"]}')
        assert response.status_code == 200
        
        # Verify deleted
        get_response = client.get(f'/api/cores/{created["core_id"]}')
        assert get_response.status_code == 404
