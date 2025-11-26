/**
 * TypeScript types for ExpOpt frontend
 */

export interface Fragment {
  fragment_id: number;
  smiles: string;
  name?: string;
  properties: Record<string, number | string>;
  attachment_points: number;
  created_at?: string;
  updated_at?: string;
}

export interface Core {
  core_id: number;
  smiles: string;
  name?: string;
  description?: string;
  properties: Record<string, number | string>;
  created_at?: string;
  updated_at?: string;
}

export interface TaskResult {
  result_id: number;
  task_id: number;
  output_smiles: string;
  score?: number;
  properties: Record<string, number | string>;
  rank?: number;
}

export interface Task {
  task_id: number;
  task_type: 'rgroup' | 'core_hopping';
  input_smiles: string;
  constraints: Record<string, number | string>;
  status: 'pending' | 'running' | 'completed' | 'failed';
  created_at?: string;
  completed_at?: string;
  error_message?: string;
  results?: TaskResult[];
}

export interface PaginatedResponse<T> {
  total: number;
  limit: number;
  offset: number;
  [key: string]: T[] | number;
}

export interface TaskSubmission {
  task_type: 'rgroup' | 'core_hopping';
  input_smiles: string;
  constraints?: Record<string, number | string>;
}
