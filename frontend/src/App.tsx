import React from 'react';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Navbar from './components/Navbar';
import Home from './pages/Home';
import TaskSubmit from './pages/TaskSubmit';
import Results from './pages/Results';
import ResultsList from './pages/ResultsList';
import Database from './pages/Database';
import './App.css';

function App() {
  return (
    <Router>
      <div className="App">
        <Navbar />
        <main className="main-content">
          <Routes>
            <Route path="/" element={<Home />} />
            <Route path="/submit" element={<TaskSubmit />} />
            <Route path="/results" element={<ResultsList />} />
            <Route path="/results/:taskId" element={<Results />} />
            <Route path="/database" element={<Database />} />
          </Routes>
        </main>
      </div>
    </Router>
  );
}

export default App;
