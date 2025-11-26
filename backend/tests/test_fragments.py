"""
Tests for Fragment API endpoints.
"""

import json
import pytest


class TestFragmentsAPI:
    """Test cases for fragments API."""
    
    def test_get_fragments_empty(self, client):
        """Test getting fragments when database is empty."""
        response = client.get('/api/fragments')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['fragments'] == []
        assert data['total'] == 0
    
    def test_create_fragment(self, client):
        """Test creating a new fragment."""
        fragment_data = {
            'smiles': 'C1=CC=CC=C1',
            'name': 'Benzene',
            'attachment_points': 1,
            'properties': {'logP': 2.0, 'MW': 78.11}
        }
        response = client.post(
            '/api/fragments',
            data=json.dumps(fragment_data),
            content_type='application/json'
        )
        assert response.status_code == 201
        data = json.loads(response.data)
        assert data['smiles'] == 'C1=CC=CC=C1'
        assert data['name'] == 'Benzene'
        assert data['properties']['logP'] == 2.0
    
    def test_create_fragment_missing_smiles(self, client):
        """Test creating fragment without SMILES returns error."""
        response = client.post(
            '/api/fragments',
            data=json.dumps({'name': 'Test'}),
            content_type='application/json'
        )
        assert response.status_code == 400
    
    def test_get_fragment_by_id(self, client):
        """Test getting a specific fragment."""
        # Create fragment first
        fragment_data = {'smiles': 'CCO', 'name': 'Ethanol'}
        create_response = client.post(
            '/api/fragments',
            data=json.dumps(fragment_data),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Get by ID
        response = client.get(f'/api/fragments/{created["fragment_id"]}')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['smiles'] == 'CCO'
    
    def test_update_fragment(self, client):
        """Test updating a fragment."""
        # Create fragment first
        create_response = client.post(
            '/api/fragments',
            data=json.dumps({'smiles': 'C', 'name': 'Methane'}),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Update
        response = client.put(
            f'/api/fragments/{created["fragment_id"]}',
            data=json.dumps({'name': 'Updated Methane'}),
            content_type='application/json'
        )
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['name'] == 'Updated Methane'
    
    def test_delete_fragment(self, client):
        """Test deleting a fragment."""
        # Create fragment first
        create_response = client.post(
            '/api/fragments',
            data=json.dumps({'smiles': 'N'}),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Delete
        response = client.delete(f'/api/fragments/{created["fragment_id"]}')
        assert response.status_code == 200
        
        # Verify deleted
        get_response = client.get(f'/api/fragments/{created["fragment_id"]}')
        assert get_response.status_code == 404
