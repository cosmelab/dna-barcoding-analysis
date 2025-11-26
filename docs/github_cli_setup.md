# GitHub CLI Authentication Setup

This guide is for **instructors and developers** who need to interact with GitHub APIs, manage packages, or use the `gh` CLI tool.

**Students do NOT need this** - they only need `git clone` and `docker pull`.

---

## Why You Need a GitHub Token

GitHub requires authentication for:
- Publishing packages to GitHub Container Registry (ghcr.io)
- Using GitHub CLI (`gh`) with packages
- Accessing private repositories via API
- Managing repository settings via CLI

---

## Quick Setup (Recommended)

### 1. Create a Personal Access Token

Go to: https://github.com/settings/tokens/new

**Token settings:**
- **Name**: `gh-cli-packages` (or any descriptive name)
- **Expiration**: 90 days (or your preference)
- **Scopes** (check these):
  - ✅ `repo` (Full control of private repositories)
  - ✅ `read:packages` (Download packages)
  - ✅ `write:packages` (Upload packages)
  - ✅ `delete:packages` (Delete packages - optional)

Click **"Generate token"** and **copy it immediately** (you won't see it again!)

### 2. Store Token Securely

**⚠️ NEVER store tokens in plain text in shell profiles!**

**Secure method (uses separate file with restricted permissions):**

```bash
# Create secure token file
touch ~/.github_token
chmod 600 ~/.github_token  # Only you can read it

# Add your token to the file
echo "export GITHUB_TOKEN=YOUR_TOKEN_HERE" > ~/.github_token

# Source it in ~/.zshrc (for zsh users)
echo 'source ~/.github_token 2>/dev/null' >> ~/.zshrc

# Source it in ~/.bashrc (for bash users)
echo 'source ~/.github_token 2>/dev/null' >> ~/.bashrc

# Reload your shell
source ~/.zshrc  # or source ~/.bashrc
```

**Replace `YOUR_TOKEN_HERE` with your actual token!**

### 3. Verify It Works

```bash
# Check if token is set
echo $GITHUB_TOKEN
# Should print your token (don't share this!)

# Test GitHub API access
gh api /user --jq '.login'
# Should print your GitHub username
```

---

## Alternative: Use `gh auth login` (Keychain Method)

**More secure** - stores token in system keychain instead of files:

```bash
gh auth login
```

Then:
1. Choose: **GitHub.com**
2. Choose: **HTTPS** protocol
3. Choose: **Paste an authentication token**
4. Paste your token
5. Press Enter

**Verification:**
```bash
gh auth status
# Should show: Logged in to github.com as YOUR_USERNAME
```

---

## Using the Token

### With GitHub CLI

```bash
# List your packages
gh api /users/YOUR_USERNAME/packages/container/PACKAGE_NAME

# Check repository details
gh repo view

# Create a release
gh release create v1.0.0 --title "Version 1.0.0" --notes "Release notes"
```

### With Docker

```bash
# Login to GitHub Container Registry
echo $GITHUB_TOKEN | docker login ghcr.io -u YOUR_USERNAME --password-stdin

# Pull from GitHub Packages
docker pull ghcr.io/YOUR_USERNAME/PACKAGE_NAME:latest

# Push to GitHub Packages
docker push ghcr.io/YOUR_USERNAME/PACKAGE_NAME:latest
```

### With curl (Direct API calls)

```bash
# Get user info
curl -H "Authorization: token $GITHUB_TOKEN" https://api.github.com/user

# List packages
curl -H "Authorization: token $GITHUB_TOKEN" \
  https://api.github.com/users/YOUR_USERNAME/packages
```

---

## Security Best Practices

### ✅ DO:
- Store tokens in secure files with `chmod 600`
- Use system keychain when possible (`gh auth login`)
- Set expiration dates on tokens
- Revoke tokens you're not using
- Use separate tokens for different purposes

### ❌ DON'T:
- Share tokens in chat, email, or screenshots
- Commit tokens to git repositories
- Store tokens in plain text in dotfiles (.zshrc, .bashrc)
- Use tokens with broader scopes than needed
- Keep expired or unused tokens active

---

## Troubleshooting

### Token not working

```bash
# Check if token is set
echo $GITHUB_TOKEN

# If empty, reload shell config
source ~/.zshrc  # or source ~/.bashrc

# Test API access
gh api /user
```

### Permission denied

Your token might not have the right scopes. Create a new token with the required scopes listed above.

### Token exposed accidentally

1. **Immediately revoke it**: https://github.com/settings/tokens
2. Create a new token
3. Update `~/.github_token` with the new token
4. Reload shell: `source ~/.zshrc`

---

## For This Project

This DNA barcoding analysis project uses GitHub Actions to automatically build and push Docker containers to both:

- **Docker Hub**: `docker.io/cosmelab/dna-barcoding-analysis`
- **GitHub Packages**: `ghcr.io/cosmelab/dna-barcoding-analysis`

**As an instructor, you need this token to:**
- Publish package updates
- Manage releases
- Update repository settings
- Troubleshoot CI/CD issues

**Students do NOT need this** - they just run:
```bash
docker pull cosmelab/dna-barcoding-analysis:latest
```

---

## Related Documentation

- [GitHub Personal Access Tokens Documentation](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens)
- [GitHub CLI Authentication](https://cli.github.com/manual/gh_auth_login)
- [GitHub Packages Documentation](https://docs.github.com/en/packages)

---

**Last Updated**: November 25, 2025
