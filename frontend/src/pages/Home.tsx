/**
 * Home page component
 */

import React from 'react';
import { Link } from 'react-router-dom';
import './Home.css';

const Home: React.FC = () => {
  return (
    <div className="home-page">
      <header className="hero">
        <h1>ExpOpt - Molecular Optimization Platform</h1>
        <p>A powerful tool for RGroup replacement and Core Hopping analysis in drug discovery</p>
      </header>

      <section className="features">
        <div className="feature-card">
          <h3>🔬 RGroup Replacement</h3>
          <p>
            Explore chemical space by replacing functional groups in your molecules
            with optimized alternatives from our curated database.
          </p>
          <Link to="/submit?type=rgroup" className="btn btn-primary">
            Start RGroup Analysis
          </Link>
        </div>

        <div className="feature-card">
          <h3>🧪 Core Hopping</h3>
          <p>
            Discover new scaffolds while maintaining desired pharmacological properties.
            Find novel structures with similar activity profiles.
          </p>
          <Link to="/submit?type=core_hopping" className="btn btn-primary">
            Start Core Hopping
          </Link>
        </div>

        <div className="feature-card">
          <h3>📊 Database Management</h3>
          <p>
            Browse and manage the fragment and core databases. Add new entries,
            search existing structures, and export data.
          </p>
          <Link to="/database" className="btn btn-secondary">
            Manage Database
          </Link>
        </div>
      </section>

      <section className="getting-started">
        <h2>Getting Started</h2>
        <ol>
          <li>Enter your molecule's SMILES notation</li>
          <li>Select the analysis type (RGroup or Core Hopping)</li>
          <li>Configure search constraints</li>
          <li>Review and analyze results</li>
        </ol>
      </section>
    </div>
  );
};

export default Home;
