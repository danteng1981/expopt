/**
 * API service for communicating with the ExpOpt backend
 */

import axios from 'axios';
import { Task, TaskSubmission, Fragment, Core } from '../types';

const API_BASE_URL = process.env.REACT_APP_API_URL || 'http://localhost:5000/api';

const api = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
});

// Task APIs
export const taskService = {
  getAll: async (params?: { task_type?: string; status?: string; limit?: number; offset?: number }) => {
    const response = await api.get('/tasks', { params });
    return response.data;
  },

  getById: async (taskId: number) => {
    const response = await api.get(`/tasks/${taskId}`);
    return response.data as Task;
  },

  create: async (submission: TaskSubmission) => {
    const response = await api.post('/tasks', submission);
    return response.data as Task;
  },

  delete: async (taskId: number) => {
    const response = await api.delete(`/tasks/${taskId}`);
    return response.data;
  },

  getResults: async (taskId: number) => {
    const response = await api.get(`/tasks/${taskId}/results`);
    return response.data;
  },
};

// Fragment APIs
export const fragmentService = {
  getAll: async (params?: { attachment_points?: number; limit?: number; offset?: number }) => {
    const response = await api.get('/fragments', { params });
    return response.data;
  },

  getById: async (fragmentId: number) => {
    const response = await api.get(`/fragments/${fragmentId}`);
    return response.data as Fragment;
  },

  create: async (fragment: Partial<Fragment>) => {
    const response = await api.post('/fragments', fragment);
    return response.data as Fragment;
  },

  update: async (fragmentId: number, fragment: Partial<Fragment>) => {
    const response = await api.put(`/fragments/${fragmentId}`, fragment);
    return response.data as Fragment;
  },

  delete: async (fragmentId: number) => {
    const response = await api.delete(`/fragments/${fragmentId}`);
    return response.data;
  },

  search: async (smiles: string, limit?: number) => {
    const response = await api.get('/fragments/search', { params: { smiles, limit } });
    return response.data;
  },
};

// Core APIs
export const coreService = {
  getAll: async (params?: { limit?: number; offset?: number }) => {
    const response = await api.get('/cores', { params });
    return response.data;
  },

  getById: async (coreId: number) => {
    const response = await api.get(`/cores/${coreId}`);
    return response.data as Core;
  },

  create: async (core: Partial<Core>) => {
    const response = await api.post('/cores', core);
    return response.data as Core;
  },

  update: async (coreId: number, core: Partial<Core>) => {
    const response = await api.put(`/cores/${coreId}`, core);
    return response.data as Core;
  },

  delete: async (coreId: number) => {
    const response = await api.delete(`/cores/${coreId}`);
    return response.data;
  },

  search: async (params: { smiles?: string; name?: string; limit?: number }) => {
    const response = await api.get('/cores/search', { params });
    return response.data;
  },
};

export default api;
