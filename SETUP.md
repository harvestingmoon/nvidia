# Setup Guide

## Environment Configuration

This project uses environment variables to securely store API keys and configuration.

### Quick Setup

1. **Copy the environment template:**
   ```bash
   cp .env.example .env
   ```

2. **Edit `.env` and add your NVIDIA API key:**
   ```bash
   # Open in your editor
   nano .env
   # or
   code .env
   ```

3. **Get your API key:**
   - Visit [NVIDIA Build](https://build.nvidia.com/)
   - Sign in with your NVIDIA account
   - Navigate to any model and click "Get API Key"
   - Copy your API key (starts with `nvapi-`)

4. **Update your `.env` file:**
   ```env
   NVIDIA_API_KEY=nvapi-YOUR_ACTUAL_KEY_HERE
   NGC_CLI_API_KEY=nvapi-YOUR_ACTUAL_KEY_HERE
   ```

5. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

6. **Run the application:**
   ```bash
   bash scripts/launch_workflow.sh
   ```

### Security Notes

- ✅ The `.env` file is in `.gitignore` - it will NOT be committed to git
- ✅ Never hardcode API keys in your source code
- ✅ Never commit the `.env` file to version control
- ✅ Use `.env.example` as a template (without real keys)

### Environment Variables

| Variable | Description | Required |
|----------|-------------|----------|
| `NVIDIA_API_KEY` | Your NVIDIA NGC API key | Yes |
| `NGC_CLI_API_KEY` | Alternative name for API key | Optional |
| `NIM_ENDPOINT` | NVIDIA NIM endpoint URL | No (has default) |
| `STREAMLIT_SERVER_PORT` | Port for Streamlit app | No (default: 8501) |

### Troubleshooting

**"No API key configured" warning:**
- Make sure `.env` file exists in the project root
- Check that `NVIDIA_API_KEY` is set in `.env`
- Verify the key starts with `nvapi-`

**Import errors for `dotenv`:**
```bash
pip install python-dotenv
```

**API key not loading:**
- Restart your application after editing `.env`
- Check file permissions on `.env`
- Make sure `.env` is in the project root directory

### For Public Repositories

When sharing your code publicly:

1. ✅ **DO** commit `.env.example` (template without keys)
2. ✅ **DO** include `.env` in `.gitignore`
3. ✅ **DO** document environment variables in README
4. ❌ **DON'T** commit `.env` with real keys
5. ❌ **DON'T** hardcode API keys in source files
6. ❌ **DON'T** include keys in documentation examples

### Verifying Your Setup

Check if environment variables are loaded:

```bash
python -c "from dotenv import load_dotenv; import os; load_dotenv(); print('API Key loaded:', bool(os.getenv('NVIDIA_API_KEY')))"
```

Should output: `API Key loaded: True`
