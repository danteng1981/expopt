/**
 * Navigation bar component
 */

import React from 'react';
import { Link, useLocation } from 'react-router-dom';
import './Navbar.css';

const Navbar: React.FC = () => {
  const location = useLocation();

  return (
    <nav className="navbar">
      <div className="navbar-brand">
        <Link to="/">ExpOpt</Link>
      </div>
      <ul className="navbar-nav">
        <li className={location.pathname === '/' ? 'active' : ''}>
          <Link to="/">Home</Link>
        </li>
        <li className={location.pathname === '/submit' ? 'active' : ''}>
          <Link to="/submit">Submit Task</Link>
        </li>
        <li className={location.pathname === '/results' ? 'active' : ''}>
          <Link to="/results">Results</Link>
        </li>
        <li className={location.pathname === '/database' ? 'active' : ''}>
          <Link to="/database">Database</Link>
        </li>
      </ul>
    </nav>
  );
};

export default Navbar;
