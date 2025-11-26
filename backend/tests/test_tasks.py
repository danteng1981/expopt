"""
Tests for Task API endpoints.
"""

import json
import pytest


class TestTasksAPI:
    """Test cases for tasks API."""
    
    def test_get_tasks_empty(self, client):
        """Test getting tasks when database is empty."""
        response = client.get('/api/tasks')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['tasks'] == []
        assert data['total'] == 0
    
    def test_create_rgroup_task(self, client):
        """Test creating an RGroup replacement task."""
        task_data = {
            'task_type': 'rgroup',
            'input_smiles': 'C1=CC=CC=C1',
            'constraints': {'max_results': 10}
        }
        response = client.post(
            '/api/tasks',
            data=json.dumps(task_data),
            content_type='application/json'
        )
        assert response.status_code == 201
        data = json.loads(response.data)
        assert data['task_type'] == 'rgroup'
        assert data['status'] == 'completed'
        assert len(data['results']) > 0
    
    def test_create_core_hopping_task(self, client):
        """Test creating a CoreHopping task."""
        task_data = {
            'task_type': 'core_hopping',
            'input_smiles': 'c1ccccc1'
        }
        response = client.post(
            '/api/tasks',
            data=json.dumps(task_data),
            content_type='application/json'
        )
        assert response.status_code == 201
        data = json.loads(response.data)
        assert data['task_type'] == 'core_hopping'
        assert data['status'] == 'completed'
    
    def test_create_task_invalid_type(self, client):
        """Test creating task with invalid type returns error."""
        task_data = {
            'task_type': 'invalid',
            'input_smiles': 'C'
        }
        response = client.post(
            '/api/tasks',
            data=json.dumps(task_data),
            content_type='application/json'
        )
        assert response.status_code == 400
    
    def test_create_task_missing_smiles(self, client):
        """Test creating task without SMILES returns error."""
        response = client.post(
            '/api/tasks',
            data=json.dumps({'task_type': 'rgroup'}),
            content_type='application/json'
        )
        assert response.status_code == 400
    
    def test_get_task_by_id(self, client):
        """Test getting a specific task with results."""
        # Create task first
        task_data = {
            'task_type': 'rgroup',
            'input_smiles': 'CCO'
        }
        create_response = client.post(
            '/api/tasks',
            data=json.dumps(task_data),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Get by ID
        response = client.get(f'/api/tasks/{created["task_id"]}')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['task_id'] == created['task_id']
        assert 'results' in data
    
    def test_get_task_results(self, client):
        """Test getting task results."""
        # Create task first
        task_data = {
            'task_type': 'rgroup',
            'input_smiles': 'C'
        }
        create_response = client.post(
            '/api/tasks',
            data=json.dumps(task_data),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Get results
        response = client.get(f'/api/tasks/{created["task_id"]}/results')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['task_id'] == created['task_id']
        assert 'results' in data
    
    def test_delete_task(self, client):
        """Test deleting a task."""
        # Create task first
        task_data = {
            'task_type': 'rgroup',
            'input_smiles': 'N'
        }
        create_response = client.post(
            '/api/tasks',
            data=json.dumps(task_data),
            content_type='application/json'
        )
        created = json.loads(create_response.data)
        
        # Delete
        response = client.delete(f'/api/tasks/{created["task_id"]}')
        assert response.status_code == 200
        
        # Verify deleted
        get_response = client.get(f'/api/tasks/{created["task_id"]}')
        assert get_response.status_code == 404
