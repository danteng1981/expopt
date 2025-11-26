/**
 * Results view page component
 */

import React, { useState, useEffect } from 'react';
import { useParams, Link } from 'react-router-dom';
import { taskService } from '../services/api';
import { Task } from '../types';
import MoleculeViewer from '../components/MoleculeViewer';
import './Results.css';

const Results: React.FC = () => {
  const { taskId } = useParams<{ taskId: string }>();
  const [task, setTask] = useState<Task | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const fetchTask = async () => {
      if (!taskId) return;
      
      try {
        const data = await taskService.getById(parseInt(taskId));
        setTask(data);
      } catch (err) {
        setError('Failed to load task results');
      } finally {
        setLoading(false);
      }
    };

    fetchTask();
  }, [taskId]);

  if (loading) {
    return (
      <div className="results-page">
        <div className="loading">Loading results...</div>
      </div>
    );
  }

  if (error || !task) {
    return (
      <div className="results-page">
        <div className="error">{error || 'Task not found'}</div>
        <Link to="/results" className="back-link">← Back to all tasks</Link>
      </div>
    );
  }

  return (
    <div className="results-page">
      <div className="results-header">
        <div className="header-content">
          <h1>Task Results</h1>
          <Link to="/results" className="back-link">← Back to all tasks</Link>
        </div>
        
        <div className="task-info">
          <div className="info-item">
            <span className="label">Task ID:</span>
            <span className="value">{task.task_id}</span>
          </div>
          <div className="info-item">
            <span className="label">Type:</span>
            <span className="value task-type">{task.task_type === 'rgroup' ? 'RGroup Replacement' : 'Core Hopping'}</span>
          </div>
          <div className="info-item">
            <span className="label">Status:</span>
            <span className={`value status status-${task.status}`}>{task.status}</span>
          </div>
          <div className="info-item">
            <span className="label">Created:</span>
            <span className="value">{task.created_at ? new Date(task.created_at).toLocaleString() : 'N/A'}</span>
          </div>
        </div>
      </div>

      <div className="input-molecule">
        <h2>Input Molecule</h2>
        <div className="molecule-display">
          <MoleculeViewer smiles={task.input_smiles} width={250} height={180} />
          <div className="smiles-text">
            <code>{task.input_smiles}</code>
          </div>
        </div>
      </div>

      {task.status === 'completed' && task.results && task.results.length > 0 ? (
        <div className="results-section">
          <h2>Results ({task.results.length})</h2>
          <div className="results-table-container">
            <table className="results-table">
              <thead>
                <tr>
                  <th>Rank</th>
                  <th>Structure</th>
                  <th>SMILES</th>
                  <th>Score</th>
                  <th>Properties</th>
                </tr>
              </thead>
              <tbody>
                {task.results.map((result) => (
                  <tr key={result.result_id}>
                    <td className="rank-cell">{result.rank}</td>
                    <td className="structure-cell">
                      <MoleculeViewer smiles={result.output_smiles} width={150} height={100} />
                    </td>
                    <td className="smiles-cell">
                      <code>{result.output_smiles}</code>
                    </td>
                    <td className="score-cell">
                      {result.score ? (
                        <span className="score-badge">{result.score.toFixed(2)}</span>
                      ) : 'N/A'}
                    </td>
                    <td className="properties-cell">
                      {Object.entries(result.properties).map(([key, value]) => (
                        <div key={key} className="property-item">
                          <span className="prop-key">{key}:</span>
                          <span className="prop-value">{typeof value === 'number' ? value.toFixed(2) : value}</span>
                        </div>
                      ))}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      ) : task.status === 'failed' ? (
        <div className="error-section">
          <h2>Error</h2>
          <p>{task.error_message || 'An unknown error occurred'}</p>
        </div>
      ) : task.status === 'pending' || task.status === 'running' ? (
        <div className="pending-section">
          <h2>Processing</h2>
          <p>Your task is being processed. Please check back later.</p>
        </div>
      ) : (
        <div className="no-results">
          <h2>No Results</h2>
          <p>No results were found for this task.</p>
        </div>
      )}
    </div>
  );
};

export default Results;
