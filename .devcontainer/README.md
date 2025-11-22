# Dev Container Configuration

This directory contains configuration for **VS Code Dev Containers**.

## What is This?

Dev Containers let you open this project **inside** the Docker container, so:
- Your VS Code terminal runs directly in the container
- No need to type `docker run` commands
- All tools (BioPython, MAFFT, IQ-TREE) are immediately available
- One-click setup

## Who Should Use This?

**For Beginners:** 
- Stick with the regular workflow (running `./tutorial.sh` and `./run-analysis.sh`)
- You don't need Dev Containers to complete the assignment

**For Advanced Users:**
- If you want to explore the container internals
- If you want to modify Python scripts and test them immediately
- If you're comfortable with Docker concepts

## How to Use

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

- Python 3.10 environment
- BioPython for sequence analysis
- MAFFT for alignment
- IQ-TREE for phylogenetics
- All dependencies pre-installed

## Troubleshooting

**Container won't start:**
- Make sure Docker Desktop is running
- Try "Rebuild Container" from Command Palette (Cmd/Ctrl+Shift+P)

**Extensions not working:**
- They install automatically on first launch
- If missing, reload VS Code window

**Want to go back to normal mode:**
- Click green button in bottom-left
- Select "Reopen Folder Locally"

## More Info

- [VS Code Dev Containers Documentation](https://code.visualstudio.com/docs/devcontainers/containers)
- [ENTM201L CLI Tools Setup](https://cosmelab.github.io/entm201l-fall2025/setup/cli-tools.html)
