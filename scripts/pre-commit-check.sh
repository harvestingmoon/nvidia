#!/bin/bash
# Pre-commit hook to prevent committing sensitive data
# To install: cp scripts/pre-commit-check.sh .git/hooks/pre-commit && chmod +x .git/hooks/pre-commit

# Colors for output
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

echo "🔍 Running pre-commit security checks..."

# Check if .env file is being committed
if git diff --cached --name-only | grep -q "^\.env$"; then
    echo -e "${RED}❌ ERROR: Attempting to commit .env file!${NC}"
    echo "   The .env file contains sensitive API keys and should never be committed."
    echo "   Please unstage it with: git reset HEAD .env"
    exit 1
fi

# Check for hardcoded API keys in staged files
if git diff --cached | grep -i "nvapi-[a-zA-Z0-9_-]\{30,\}"; then
    echo -e "${RED}❌ ERROR: Found potential API key in staged changes!${NC}"
    echo "   Please remove hardcoded API keys and use environment variables instead."
    echo "   Example: os.getenv('NVIDIA_API_KEY')"
    exit 1
fi

# Check for other common secrets
if git diff --cached | grep -iE "(password|secret|token).*=.*['\"][^'\"]{8,}"; then
    echo -e "${YELLOW}⚠️  WARNING: Found potential hardcoded secret!${NC}"
    echo "   Please review your changes for sensitive data."
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo -e "${GREEN}✅ Security checks passed${NC}"
exit 0
