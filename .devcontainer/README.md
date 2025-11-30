# Dev Container Configuration

This directory contains configuration for **VS Code Dev Containers** and **GitHub Codespaces**.

## What is This?

Dev Containers and Codespaces let you open this project **inside** the Docker container, so:
- Your VS Code terminal runs directly in the container
- No need to type `docker run` commands
- All tools (BioPython, MAFFT, IQ-TREE) are immediately available
- One-click setup

## Who Should Use This?

**For Beginners:**
- Stick with the regular workflow (running `./tutorial.sh` and `./run-analysis.sh`)
- You don't need Dev Containers or Codespaces to complete the assignment
- **GitHub Codespaces** is great if you don't want to install Docker locally

**For Advanced Users:**
- If you want to explore the container internals
- If you want to modify Python scripts and test them immediately
- If you're comfortable with Docker concepts

## Option 1: GitHub Codespaces (Easiest)

**Best for:** Students who don't want to install Docker, or want to work from any computer.

1. **Open in Codespaces:**
   - Go to your GitHub repository page
   - Click the green "Code" button
   - Select "Codespaces" tab
   - Click "Create codespace on main"
   - Wait ~1-2 minutes for setup

2. **Start working:**
   - Terminal is ready in the browser-based VS Code
   - Run: `./tutorial.sh` (no Docker installation needed!)
   - Run: `./run-analysis.sh` for the main analysis

3. **View HTML Reports:**

   You have **three options** to view the HTML reports generated in `results/`:

   **Option A: Live Server (Easiest)**
   - Right-click any `.html` file in the `results/` folder
   - Select "Open with Live Server"
   - Report opens in a new browser tab automatically

   **Option B: Python HTTP Server**
   ```bash
   # Start a web server in the results directory
   python -m http.server 8000 -d results/
   ```
   - Click the notification to open port 8000 in browser
   - Or go to "Ports" tab and click the globe icon next to port 8000
   - Navigate to the HTML file you want to view

   **Option C: Download and View Locally**
   - Right-click any file in `results/` folder
   - Select "Download"
   - Open the downloaded HTML file in your browser

4. **Stop/Resume Codespace:**
   - Codespaces auto-stop after 30 min of inactivity (saves your free hours)
   - Resume anytime from github.com (your work is saved)
   - Free tier: 60 hours/month for students

5. **Submit your work:**
   ```bash
   git add answers.json results/
   git commit -m "Complete assignment"
   git push
   ```

## Option 2: Local Dev Containers

**Best for:** Advanced users who want full control and faster performance.

1. **Install Dev Containers extension** (you should already have it from ENTM201L):
   - Extension ID: `ms-vscode-remote.remote-containers`
   - Install from VS Code Extensions marketplace

2. **Open this project in VS Code:**
   ```bash
   code .
   ```

3. **Reopen in Container:**
   - Click the green button in bottom-left corner of VS Code
   - Select "Reopen in Container"
   - Wait for container to start (~30 seconds first time)

4. **Start working:**
   - Terminal is now inside the container
   - Just run: `./tutorial.sh` (no `docker run` needed!)

## What's Included?

Both Codespaces and local Dev Containers provide:

- Python 3.10 environment
- BioPython for sequence analysis
- MAFFT for alignment
- IQ-TREE for phylogenetics
- All dependencies pre-installed
- VS Code extensions:
  - Python language support
  - Live Server for viewing HTML reports
  - Markdown tools
  - Docker support (local only)

## Troubleshooting

### GitHub Codespaces Issues

**Can't see HTML reports:**
- Use Live Server: Right-click .html file → "Open with Live Server"
- Or use Python server: `python -m http.server 8000 -d results/`
- Check the "Ports" tab to see forwarded ports
- If all else fails, download the HTML file and open locally

**Codespace is slow:**
- Free tier uses 2-core machines (slower than local)
- Consider upgrading to 4-core if you have GitHub Pro/Student benefits
- Or switch to local Dev Container for better performance

**Running out of free hours:**
- Free tier: 60 hours/month
- Stop your Codespace when not using it (auto-stops after 30 min)
- Use local Dev Containers if you need unlimited time

**Can't push to GitHub:**
- Make sure you're authenticated with GitHub
- Use the integrated Source Control tab in VS Code
- Or run: `gh auth login` if using GitHub CLI

### Local Dev Container Issues

**Container won't start:**
- Make sure Docker Desktop is running
- Try "Rebuild Container" from Command Palette (Cmd/Ctrl+Shift+P)

**Extensions not working:**
- They install automatically on first launch
- If missing, reload VS Code window

**Want to go back to normal mode:**
- Click green button in bottom-left
- Select "Reopen Folder Locally"

## Regular Workflow vs. Codespaces

**Important:** The regular workflow scripts (`./tutorial.sh` and `./run-analysis.sh`) work the same in both environments:

- **On your local computer:** Scripts use `docker run` automatically
- **In Codespaces/Dev Containers:** Scripts run directly (already inside container)
- No special `-cs` suffix scripts are needed - the same scripts work everywhere!

The pipeline automatically detects the environment and adjusts accordingly.

## More Info

- [GitHub Codespaces Documentation](https://docs.github.com/en/codespaces)
- [VS Code Dev Containers Documentation](https://code.visualstudio.com/docs/devcontainers/containers)
- [ENTM201L CLI Tools Setup](https://cosmelab.github.io/entm201l-fall2025/setup/cli-tools.html)
