"""
Fragment API routes for managing molecular fragments.
"""

from flask import Blueprint, request, jsonify
from app import db
from app.models.fragment import Fragment

fragments_bp = Blueprint('fragments', __name__)


@fragments_bp.route('', methods=['GET'])
def get_fragments():
    """Get all fragments with optional filtering."""
    # Query parameters for filtering
    attachment_points = request.args.get('attachment_points', type=int)
    limit = request.args.get('limit', 100, type=int)
    offset = request.args.get('offset', 0, type=int)
    
    query = Fragment.query
    
    if attachment_points is not None:
        query = query.filter_by(attachment_points=attachment_points)
    
    fragments = query.offset(offset).limit(limit).all()
    total = query.count()
    
    return jsonify({
        'fragments': [f.to_dict() for f in fragments],
        'total': total,
        'limit': limit,
        'offset': offset
    })


@fragments_bp.route('/<int:fragment_id>', methods=['GET'])
def get_fragment(fragment_id):
    """Get a specific fragment by ID."""
    fragment = Fragment.query.get_or_404(fragment_id)
    return jsonify(fragment.to_dict())


@fragments_bp.route('', methods=['POST'])
def create_fragment():
    """Create a new fragment."""
    data = request.get_json()
    
    if not data or 'smiles' not in data:
        return jsonify({'error': 'SMILES is required'}), 400
    
    fragment = Fragment(
        smiles=data['smiles'],
        name=data.get('name'),
        attachment_points=data.get('attachment_points', 1)
    )
    
    if 'properties' in data:
        fragment.set_properties(data['properties'])
    
    db.session.add(fragment)
    db.session.commit()
    
    return jsonify(fragment.to_dict()), 201


@fragments_bp.route('/<int:fragment_id>', methods=['PUT'])
def update_fragment(fragment_id):
    """Update an existing fragment."""
    fragment = Fragment.query.get_or_404(fragment_id)
    data = request.get_json()
    
    if 'smiles' in data:
        fragment.smiles = data['smiles']
    if 'name' in data:
        fragment.name = data['name']
    if 'attachment_points' in data:
        fragment.attachment_points = data['attachment_points']
    if 'properties' in data:
        fragment.set_properties(data['properties'])
    
    db.session.commit()
    return jsonify(fragment.to_dict())


@fragments_bp.route('/<int:fragment_id>', methods=['DELETE'])
def delete_fragment(fragment_id):
    """Delete a fragment."""
    fragment = Fragment.query.get_or_404(fragment_id)
    db.session.delete(fragment)
    db.session.commit()
    return jsonify({'message': 'Fragment deleted successfully'})


@fragments_bp.route('/search', methods=['GET'])
def search_fragments():
    """Search fragments by SMILES pattern."""
    smiles_pattern = request.args.get('smiles', '')
    limit = request.args.get('limit', 100, type=int)
    
    if not smiles_pattern:
        return jsonify({'error': 'SMILES pattern is required'}), 400
    
    fragments = Fragment.query.filter(
        Fragment.smiles.contains(smiles_pattern)
    ).limit(limit).all()
    
    return jsonify({
        'fragments': [f.to_dict() for f in fragments],
        'total': len(fragments)
    })
