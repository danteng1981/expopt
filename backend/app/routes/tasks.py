"""
Task API routes for submitting and managing replacement tasks.
"""

from datetime import datetime, timezone
from flask import Blueprint, request, jsonify
from app import db
from app.models.task import Task, TaskResult

tasks_bp = Blueprint('tasks', __name__)


@tasks_bp.route('', methods=['GET'])
def get_tasks():
    """Get all tasks with optional filtering."""
    task_type = request.args.get('task_type')
    status = request.args.get('status')
    limit = request.args.get('limit', 100, type=int)
    offset = request.args.get('offset', 0, type=int)
    
    query = Task.query
    
    if task_type:
        query = query.filter_by(task_type=task_type)
    if status:
        query = query.filter_by(status=status)
    
    tasks = query.order_by(Task.created_at.desc()).offset(offset).limit(limit).all()
    total = query.count()
    
    return jsonify({
        'tasks': [t.to_dict() for t in tasks],
        'total': total,
        'limit': limit,
        'offset': offset
    })


@tasks_bp.route('/<int:task_id>', methods=['GET'])
def get_task(task_id):
    """Get a specific task by ID with results."""
    task = Task.query.get_or_404(task_id)
    return jsonify(task.to_dict(include_results=True))


@tasks_bp.route('', methods=['POST'])
def create_task():
    """Submit a new replacement task."""
    data = request.get_json()
    
    if not data:
        return jsonify({'error': 'Request body is required'}), 400
    
    if 'task_type' not in data:
        return jsonify({'error': 'task_type is required'}), 400
    
    if data['task_type'] not in ['rgroup', 'core_hopping']:
        return jsonify({'error': 'task_type must be "rgroup" or "core_hopping"'}), 400
    
    if 'input_smiles' not in data:
        return jsonify({'error': 'input_smiles is required'}), 400
    
    task = Task(
        task_type=data['task_type'],
        input_smiles=data['input_smiles'],
        status='pending'
    )
    
    if 'constraints' in data:
        task.set_constraints(data['constraints'])
    
    db.session.add(task)
    db.session.commit()
    
    # Simulate task processing (in real implementation, this would be async)
    _process_task(task)
    
    return jsonify(task.to_dict(include_results=True)), 201


def _process_task(task):
    """Process a task and generate mock results."""
    try:
        task.status = 'running'
        db.session.commit()
        
        # Generate mock results based on task type
        if task.task_type == 'rgroup':
            results = _generate_rgroup_results(task)
        else:
            results = _generate_core_hopping_results(task)
        
        for result_data in results:
            result = TaskResult(
                task_id=task.task_id,
                output_smiles=result_data['smiles'],
                score=result_data.get('score'),
                rank=result_data.get('rank')
            )
            if 'properties' in result_data:
                result.set_properties(result_data['properties'])
            db.session.add(result)
        
        task.status = 'completed'
        task.completed_at = datetime.now(timezone.utc)
        db.session.commit()
        
    except Exception as e:
        task.status = 'failed'
        task.error_message = str(e)
        db.session.commit()


def _generate_rgroup_results(task):
    """Generate mock RGroup replacement results."""
    base_smiles = task.input_smiles
    
    # Mock results - in real implementation, this would use RDKit
    return [
        {
            'smiles': f'{base_smiles}-CH3',
            'score': 0.95,
            'rank': 1,
            'properties': {'logP': 2.1, 'MW': 150.2}
        },
        {
            'smiles': f'{base_smiles}-CF3',
            'score': 0.88,
            'rank': 2,
            'properties': {'logP': 2.8, 'MW': 186.2}
        },
        {
            'smiles': f'{base_smiles}-OCH3',
            'score': 0.82,
            'rank': 3,
            'properties': {'logP': 1.5, 'MW': 166.2}
        }
    ]


def _generate_core_hopping_results(task):
    """Generate mock CoreHopping results."""
    return [
        {
            'smiles': 'c1ccc2[nH]ccc2c1',
            'score': 0.91,
            'rank': 1,
            'properties': {'similarity': 0.85, 'MW': 117.1}
        },
        {
            'smiles': 'c1ccc2occc2c1',
            'score': 0.85,
            'rank': 2,
            'properties': {'similarity': 0.78, 'MW': 118.1}
        }
    ]


@tasks_bp.route('/<int:task_id>', methods=['DELETE'])
def delete_task(task_id):
    """Delete a task and its results."""
    task = Task.query.get_or_404(task_id)
    db.session.delete(task)
    db.session.commit()
    return jsonify({'message': 'Task deleted successfully'})


@tasks_bp.route('/<int:task_id>/results', methods=['GET'])
def get_task_results(task_id):
    """Get results for a specific task."""
    task = Task.query.get_or_404(task_id)
    results = TaskResult.query.filter_by(task_id=task_id).order_by(TaskResult.rank).all()
    
    return jsonify({
        'task_id': task_id,
        'status': task.status,
        'results': [r.to_dict() for r in results]
    })
