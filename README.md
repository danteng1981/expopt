# ExpOpt

Molecular Optimization Platform for medicine and chemistry.

## Overview

ExpOpt is a web application for RGroup replacement and Core Hopping analysis in drug discovery. It provides:

- **Database Support**: SQLAlchemy-based storage for fragments, cores, and task results
- **REST API**: Flask backend with endpoints for task submission and result retrieval
- **Frontend Interface**: React-based UI for submitting tasks and visualizing results

## Project Structure

```
expopt/
├── backend/           # Flask backend application
│   ├── app/
│   │   ├── models/    # SQLAlchemy database models
│   │   ├── routes/    # API endpoints
│   │   └── services/  # Business logic
│   ├── tests/         # Backend unit tests
│   ├── requirements.txt
│   └── run.py         # Application entry point
├── frontend/          # React frontend application
│   ├── src/
│   │   ├── components/  # Reusable UI components
│   │   ├── pages/       # Page components
│   │   ├── services/    # API client services
│   │   └── types/       # TypeScript type definitions
│   └── package.json
└── README.md
```

## Getting Started

### Backend Setup

1. Navigate to the backend directory:
   ```bash
   cd backend
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Run the Flask server:
   ```bash
   python run.py
   ```

The backend API will be available at `http://localhost:5000`.

### Frontend Setup

1. Navigate to the frontend directory:
   ```bash
   cd frontend
   ```

2. Install dependencies:
   ```bash
   npm install
   ```

3. Start the development server:
   ```bash
   npm start
   ```

The frontend will be available at `http://localhost:3000`.

## API Endpoints

### Tasks

- `GET /api/tasks` - List all tasks
- `POST /api/tasks` - Submit a new task
- `GET /api/tasks/<id>` - Get task details with results
- `DELETE /api/tasks/<id>` - Delete a task
- `GET /api/tasks/<id>/results` - Get task results

### Fragments

- `GET /api/fragments` - List all fragments
- `POST /api/fragments` - Create a fragment
- `GET /api/fragments/<id>` - Get fragment details
- `PUT /api/fragments/<id>` - Update a fragment
- `DELETE /api/fragments/<id>` - Delete a fragment
- `GET /api/fragments/search` - Search fragments

### Cores

- `GET /api/cores` - List all cores
- `POST /api/cores` - Create a core
- `GET /api/cores/<id>` - Get core details
- `PUT /api/cores/<id>` - Update a core
- `DELETE /api/cores/<id>` - Delete a core
- `GET /api/cores/search` - Search cores

## Database Schema

### Fragments Table
- `fragment_id` (primary key)
- `smiles` - SMILES structure
- `properties` - JSON-encoded physicochemical attributes
- `attachment_points` - Number of attachment points
- `name` - Optional name

### Cores Table
- `core_id` (primary key)
- `smiles` - SMILES structure/scaffold
- `name` - Optional name
- `description` - Optional description
- `properties` - JSON-encoded properties

### Tasks Table
- `task_id` (primary key)
- `task_type` - 'rgroup' or 'core_hopping'
- `input_smiles` - Input molecule SMILES
- `constraints` - JSON-encoded search constraints
- `status` - pending, running, completed, failed

### Task Results Table
- `result_id` (primary key)
- `task_id` (foreign key)
- `output_smiles` - Result molecule SMILES
- `score` - Similarity/scoring value
- `properties` - JSON-encoded properties
- `rank` - Result ranking

## Running Tests

### Backend Tests
```bash
cd backend
python -m pytest tests/ -v
```

### Frontend Tests
```bash
cd frontend
npm test
```

## License

MIT License
