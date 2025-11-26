"""
Core API routes for managing molecular cores/scaffolds.
"""

from flask import Blueprint, request, jsonify
from app import db
from app.models.core import Core

cores_bp = Blueprint('cores', __name__)


@cores_bp.route('', methods=['GET'])
def get_cores():
    """Get all cores with optional pagination."""
    limit = request.args.get('limit', 100, type=int)
    offset = request.args.get('offset', 0, type=int)
    
    query = Core.query
    cores = query.offset(offset).limit(limit).all()
    total = query.count()
    
    return jsonify({
        'cores': [c.to_dict() for c in cores],
        'total': total,
        'limit': limit,
        'offset': offset
    })


@cores_bp.route('/<int:core_id>', methods=['GET'])
def get_core(core_id):
    """Get a specific core by ID."""
    core = Core.query.get_or_404(core_id)
    return jsonify(core.to_dict())


@cores_bp.route('', methods=['POST'])
def create_core():
    """Create a new core."""
    data = request.get_json()
    
    if not data or 'smiles' not in data:
        return jsonify({'error': 'SMILES is required'}), 400
    
    core = Core(
        smiles=data['smiles'],
        name=data.get('name'),
        description=data.get('description')
    )
    
    if 'properties' in data:
        core.set_properties(data['properties'])
    
    db.session.add(core)
    db.session.commit()
    
    return jsonify(core.to_dict()), 201


@cores_bp.route('/<int:core_id>', methods=['PUT'])
def update_core(core_id):
    """Update an existing core."""
    core = Core.query.get_or_404(core_id)
    data = request.get_json()
    
    if 'smiles' in data:
        core.smiles = data['smiles']
    if 'name' in data:
        core.name = data['name']
    if 'description' in data:
        core.description = data['description']
    if 'properties' in data:
        core.set_properties(data['properties'])
    
    db.session.commit()
    return jsonify(core.to_dict())


@cores_bp.route('/<int:core_id>', methods=['DELETE'])
def delete_core(core_id):
    """Delete a core."""
    core = Core.query.get_or_404(core_id)
    db.session.delete(core)
    db.session.commit()
    return jsonify({'message': 'Core deleted successfully'})


@cores_bp.route('/search', methods=['GET'])
def search_cores():
    """Search cores by SMILES pattern or name."""
    smiles_pattern = request.args.get('smiles', '')
    name_pattern = request.args.get('name', '')
    limit = request.args.get('limit', 100, type=int)
    
    query = Core.query
    
    if smiles_pattern:
        query = query.filter(Core.smiles.contains(smiles_pattern))
    if name_pattern:
        query = query.filter(Core.name.contains(name_pattern))
    
    cores = query.limit(limit).all()
    
    return jsonify({
        'cores': [c.to_dict() for c in cores],
        'total': len(cores)
    })
