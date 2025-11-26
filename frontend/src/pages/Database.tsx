/**
 * Database management page component
 */

import React, { useState, useEffect, useCallback } from 'react';
import { fragmentService, coreService } from '../services/api';
import { Fragment, Core } from '../types';
import './Database.css';

type TabType = 'fragments' | 'cores';

const Database: React.FC = () => {
  const [activeTab, setActiveTab] = useState<TabType>('fragments');
  const [fragments, setFragments] = useState<Fragment[]>([]);
  const [cores, setCores] = useState<Core[]>([]);
  const [loading, setLoading] = useState(true);
  const [showAddForm, setShowAddForm] = useState(false);
  const [newItem, setNewItem] = useState({ smiles: '', name: '', attachment_points: 1 });

  const fetchData = useCallback(async () => {
    setLoading(true);
    try {
      if (activeTab === 'fragments') {
        const data = await fragmentService.getAll();
        setFragments(data.fragments || []);
      } else {
        const data = await coreService.getAll();
        setCores(data.cores || []);
      }
    } catch (err) {
      console.error('Failed to fetch data:', err);
    } finally {
      setLoading(false);
    }
  }, [activeTab]);

  useEffect(() => {
    fetchData();
  }, [fetchData]);

  const handleAdd = async (e: React.FormEvent) => {
    e.preventDefault();
    try {
      if (activeTab === 'fragments') {
        await fragmentService.create({
          smiles: newItem.smiles,
          name: newItem.name,
          attachment_points: newItem.attachment_points,
        });
      } else {
        await coreService.create({
          smiles: newItem.smiles,
          name: newItem.name,
        });
      }
      setNewItem({ smiles: '', name: '', attachment_points: 1 });
      setShowAddForm(false);
      fetchData();
    } catch (err) {
      alert('Failed to add item');
    }
  };

  const handleDelete = async (id: number) => {
    if (!window.confirm('Are you sure you want to delete this item?')) return;
    
    try {
      if (activeTab === 'fragments') {
        await fragmentService.delete(id);
      } else {
        await coreService.delete(id);
      }
      fetchData();
    } catch (err) {
      alert('Failed to delete item');
    }
  };

  return (
    <div className="database-page">
      <h1>Database Management</h1>

      <div className="tabs">
        <button
          className={`tab ${activeTab === 'fragments' ? 'active' : ''}`}
          onClick={() => setActiveTab('fragments')}
        >
          Fragments ({fragments.length})
        </button>
        <button
          className={`tab ${activeTab === 'cores' ? 'active' : ''}`}
          onClick={() => setActiveTab('cores')}
        >
          Cores ({cores.length})
        </button>
      </div>

      <div className="toolbar">
        <button className="add-btn" onClick={() => setShowAddForm(!showAddForm)}>
          {showAddForm ? 'Cancel' : `+ Add ${activeTab === 'fragments' ? 'Fragment' : 'Core'}`}
        </button>
      </div>

      {showAddForm && (
        <form onSubmit={handleAdd} className="add-form">
          <div className="form-row">
            <div className="form-group">
              <label>SMILES</label>
              <input
                type="text"
                value={newItem.smiles}
                onChange={(e) => setNewItem({ ...newItem, smiles: e.target.value })}
                placeholder="e.g., C1=CC=CC=C1"
                required
              />
            </div>
            <div className="form-group">
              <label>Name</label>
              <input
                type="text"
                value={newItem.name}
                onChange={(e) => setNewItem({ ...newItem, name: e.target.value })}
                placeholder="Optional name"
              />
            </div>
            {activeTab === 'fragments' && (
              <div className="form-group">
                <label>Attachment Points</label>
                <input
                  type="number"
                  value={newItem.attachment_points}
                  onChange={(e) => setNewItem({ ...newItem, attachment_points: parseInt(e.target.value) || 1 })}
                  min={1}
                />
              </div>
            )}
            <button type="submit" className="submit-btn">Add</button>
          </div>
        </form>
      )}

      {loading ? (
        <div className="loading">Loading...</div>
      ) : (
        <div className="data-table-container">
          <table className="data-table">
            <thead>
              <tr>
                <th>ID</th>
                <th>SMILES</th>
                <th>Name</th>
                {activeTab === 'fragments' && <th>Attachment Points</th>}
                <th>Properties</th>
                <th>Actions</th>
              </tr>
            </thead>
            <tbody>
              {activeTab === 'fragments' ? (
                fragments.map((fragment) => (
                  <tr key={fragment.fragment_id}>
                    <td>{fragment.fragment_id}</td>
                    <td className="smiles-cell"><code>{fragment.smiles}</code></td>
                    <td>{fragment.name || '-'}</td>
                    <td>{fragment.attachment_points}</td>
                    <td className="properties-cell">
                      {Object.entries(fragment.properties).map(([k, v]) => (
                        <span key={k} className="prop-tag">{k}: {v}</span>
                      ))}
                    </td>
                    <td>
                      <button
                        className="delete-btn"
                        onClick={() => handleDelete(fragment.fragment_id)}
                      >
                        Delete
                      </button>
                    </td>
                  </tr>
                ))
              ) : (
                cores.map((core) => (
                  <tr key={core.core_id}>
                    <td>{core.core_id}</td>
                    <td className="smiles-cell"><code>{core.smiles}</code></td>
                    <td>{core.name || '-'}</td>
                    <td className="properties-cell">
                      {Object.entries(core.properties).map(([k, v]) => (
                        <span key={k} className="prop-tag">{k}: {v}</span>
                      ))}
                    </td>
                    <td>
                      <button
                        className="delete-btn"
                        onClick={() => handleDelete(core.core_id)}
                      >
                        Delete
                      </button>
                    </td>
                  </tr>
                ))
              )}
            </tbody>
          </table>
          
          {((activeTab === 'fragments' && fragments.length === 0) ||
            (activeTab === 'cores' && cores.length === 0)) && (
            <div className="no-data">
              No {activeTab} found. Add some using the button above.
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default Database;
