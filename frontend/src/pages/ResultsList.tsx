/**
 * Results list page component
 */

import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import { taskService } from '../services/api';
import { Task } from '../types';
import './ResultsList.css';

const ResultsList: React.FC = () => {
  const [tasks, setTasks] = useState<Task[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [filterType, setFilterType] = useState<string>('');
  const [filterStatus, setFilterStatus] = useState<string>('');

  useEffect(() => {
    const fetchTasks = async () => {
      try {
        const params: Record<string, string> = {};
        if (filterType) params.task_type = filterType;
        if (filterStatus) params.status = filterStatus;
        
        const data = await taskService.getAll(params);
        setTasks(data.tasks || []);
      } catch (err) {
        setError('Failed to load tasks');
      } finally {
        setLoading(false);
      }
    };

    fetchTasks();
  }, [filterType, filterStatus]);

  const handleDelete = async (taskId: number) => {
    if (window.confirm('Are you sure you want to delete this task?')) {
      try {
        await taskService.delete(taskId);
        setTasks(tasks.filter(t => t.task_id !== taskId));
      } catch (err) {
        alert('Failed to delete task');
      }
    }
  };

  if (loading) {
    return (
      <div className="results-list-page">
        <div className="loading">Loading tasks...</div>
      </div>
    );
  }

  return (
    <div className="results-list-page">
      <div className="page-header">
        <h1>Task Results</h1>
        <Link to="/submit" className="new-task-btn">+ New Task</Link>
      </div>

      <div className="filters">
        <div className="filter-group">
          <label>Task Type:</label>
          <select value={filterType} onChange={(e) => setFilterType(e.target.value)}>
            <option value="">All</option>
            <option value="rgroup">RGroup Replacement</option>
            <option value="core_hopping">Core Hopping</option>
          </select>
        </div>
        
        <div className="filter-group">
          <label>Status:</label>
          <select value={filterStatus} onChange={(e) => setFilterStatus(e.target.value)}>
            <option value="">All</option>
            <option value="completed">Completed</option>
            <option value="pending">Pending</option>
            <option value="running">Running</option>
            <option value="failed">Failed</option>
          </select>
        </div>
      </div>

      {error && <div className="error-message">{error}</div>}

      {tasks.length === 0 ? (
        <div className="no-tasks">
          <p>No tasks found.</p>
          <Link to="/submit">Submit your first task</Link>
        </div>
      ) : (
        <div className="tasks-grid">
          {tasks.map((task) => (
            <div key={task.task_id} className="task-card">
              <div className="task-card-header">
                <span className={`task-badge ${task.task_type}`}>
                  {task.task_type === 'rgroup' ? '🔬 RGroup' : '🧪 Core Hopping'}
                </span>
                <span className={`status-badge status-${task.status}`}>
                  {task.status}
                </span>
              </div>
              
              <div className="task-card-body">
                <div className="task-id">Task #{task.task_id}</div>
                <div className="task-smiles" title={task.input_smiles}>
                  <code>{task.input_smiles.length > 30 ? `${task.input_smiles.substring(0, 30)}...` : task.input_smiles}</code>
                </div>
                <div className="task-date">
                  Created: {task.created_at ? new Date(task.created_at).toLocaleDateString() : 'N/A'}
                </div>
              </div>
              
              <div className="task-card-actions">
                <Link to={`/results/${task.task_id}`} className="view-btn">
                  View Results
                </Link>
                <button 
                  onClick={() => handleDelete(task.task_id)} 
                  className="delete-btn"
                >
                  Delete
                </button>
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
};

export default ResultsList;
