/**
 * Task submission page component
 */

import React, { useState } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import { taskService } from '../services/api';
import { TaskSubmission } from '../types';
import './TaskSubmit.css';

const TaskSubmit: React.FC = () => {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const initialType = searchParams.get('type') as 'rgroup' | 'core_hopping' || 'rgroup';

  const [taskType, setTaskType] = useState<'rgroup' | 'core_hopping'>(initialType);
  const [inputSmiles, setInputSmiles] = useState('');
  const [maxResults, setMaxResults] = useState(10);
  const [minScore, setMinScore] = useState(0.5);
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setError(null);

    if (!inputSmiles.trim()) {
      setError('Please enter a valid SMILES string');
      return;
    }

    setIsSubmitting(true);

    try {
      const submission: TaskSubmission = {
        task_type: taskType,
        input_smiles: inputSmiles.trim(),
        constraints: {
          max_results: maxResults,
          min_score: minScore,
        },
      };

      const result = await taskService.create(submission);
      navigate(`/results/${result.task_id}`);
    } catch (err: unknown) {
      const errorMessage = err instanceof Error ? err.message : 'Failed to submit task';
      setError(errorMessage);
    } finally {
      setIsSubmitting(false);
    }
  };

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (event) => {
        const content = event.target?.result as string;
        // Extract first SMILES from file (assuming one per line)
        const smiles = content.split('\n')[0].trim();
        setInputSmiles(smiles);
      };
      reader.readAsText(file);
    }
  };

  return (
    <div className="task-submit-page">
      <h1>Submit Analysis Task</h1>
      
      <form onSubmit={handleSubmit} className="task-form">
        <div className="form-section">
          <h2>Task Type</h2>
          <div className="task-type-selector">
            <label className={`task-type-option ${taskType === 'rgroup' ? 'selected' : ''}`}>
              <input
                type="radio"
                name="taskType"
                value="rgroup"
                checked={taskType === 'rgroup'}
                onChange={() => setTaskType('rgroup')}
              />
              <span className="option-content">
                <span className="option-icon">🔬</span>
                <span className="option-label">RGroup Replacement</span>
                <span className="option-desc">Replace functional groups with optimized alternatives</span>
              </span>
            </label>
            
            <label className={`task-type-option ${taskType === 'core_hopping' ? 'selected' : ''}`}>
              <input
                type="radio"
                name="taskType"
                value="core_hopping"
                checked={taskType === 'core_hopping'}
                onChange={() => setTaskType('core_hopping')}
              />
              <span className="option-content">
                <span className="option-icon">🧪</span>
                <span className="option-label">Core Hopping</span>
                <span className="option-desc">Find alternative scaffolds with similar properties</span>
              </span>
            </label>
          </div>
        </div>

        <div className="form-section">
          <h2>Input Molecule</h2>
          <div className="form-group">
            <label htmlFor="smiles">SMILES String</label>
            <textarea
              id="smiles"
              value={inputSmiles}
              onChange={(e) => setInputSmiles(e.target.value)}
              placeholder="Enter SMILES notation (e.g., C1=CC=CC=C1)"
              rows={3}
            />
          </div>
          
          <div className="form-group">
            <label htmlFor="file">Or upload a file</label>
            <input
              type="file"
              id="file"
              accept=".smi,.txt"
              onChange={handleFileUpload}
            />
          </div>
        </div>

        <div className="form-section">
          <h2>Constraints</h2>
          <div className="constraints-grid">
            <div className="form-group">
              <label htmlFor="maxResults">Max Results</label>
              <input
                type="number"
                id="maxResults"
                value={maxResults}
                onChange={(e) => setMaxResults(parseInt(e.target.value) || 10)}
                min={1}
                max={100}
              />
            </div>
            
            <div className="form-group">
              <label htmlFor="minScore">Min Score</label>
              <input
                type="number"
                id="minScore"
                value={minScore}
                onChange={(e) => setMinScore(parseFloat(e.target.value) || 0.5)}
                min={0}
                max={1}
                step={0.1}
              />
            </div>
          </div>
        </div>

        {error && <div className="error-message">{error}</div>}

        <button type="submit" className="submit-btn" disabled={isSubmitting}>
          {isSubmitting ? 'Submitting...' : 'Submit Task'}
        </button>
      </form>
    </div>
  );
};

export default TaskSubmit;
