# Schwab Login Automation Script

This Python script automates the login process for Charles Schwab's client portal using Selenium with Firefox.

## ⚠️ Important Security Notice

This script is for educational and personal automation purposes only. Please ensure you:
- Keep your credentials secure and never commit them to version control
- Use this script responsibly and in accordance with Schwab's terms of service
- Be aware that automated login scripts may violate some websites' terms of service
- Consider the security implications of storing credentials

## Prerequisites

1. **Python 3.7+** - Make sure Python is installed on your system
2. **Firefox Browser** - The script uses Firefox WebDriver
3. **Internet Connection** - Required to download WebDriver and access Schwab

## Installation

1. **Install Python dependencies:**
   ```bash
   pip3 install -r requirements.txt
   ```

2. **Set up environment variables (recommended):**
   ```bash
   cp .env.example .env
   # Edit .env file with your actual credentials
   ```

3. **Make the script executable:**
   ```bash
   chmod +x schwab_login.py
   ```

## Usage

### Method 1: Using Environment Variables (Recommended)
1. Edit the `.env` file with your Schwab credentials
2. Run the script:
   ```bash
   python3 schwab_login.py
   ```

### Method 2: Interactive Input
If no environment variables are set, the script will prompt for credentials:
```bash
python3 schwab_login.py
```

### Method 3: Headless Mode
To run without opening a visible browser window:
```python
# In schwab_login.py, change this line:
automator = SchwabLoginAutomator(headless=True)
```

## Features

- **Robust Element Detection**: Uses multiple selectors to find login elements
- **Error Handling**: Comprehensive error handling and logging
- **Screenshot Capture**: Takes screenshots at key points for debugging
- **Two-Factor Authentication**: Handles 2FA prompts when detected
- **Logging**: Detailed logging to both file and console
- **Security**: Supports environment variables for credential storage

## Configuration Options

You can modify the `SchwabLoginAutomator` class initialization in the `main()` function:

```python
automator = SchwabLoginAutomator(
    headless=False,  # Set to True for headless mode
    timeout=30       # Timeout in seconds for web elements
)
```

## Output Files

The script generates several output files:
- `schwab_login.log` - Detailed execution log
- `login_page.png` - Screenshot of the login page
- `after_login.png` - Screenshot after successful login
- `login_failed.png` - Screenshot if login fails

## Troubleshooting

### Common Issues

1. **WebDriver Issues:**
   ```bash
   # The script automatically downloads the latest GeckoDriver
   # If issues persist, try updating Firefox browser
   ```

2. **Element Not Found:**
   - Schwab may have updated their website structure
   - Check the screenshots to see the current page state
   - Update selectors in the `find_login_elements()` method if needed

3. **Login Fails:**
   - Verify your credentials are correct
   - Check if 2FA is required
   - Ensure your account is not locked

4. **Timeout Errors:**
   - Increase the timeout value in the constructor
   - Check your internet connection
   - Schwab's servers might be slow

### Debugging

1. **Enable Verbose Logging:**
   ```python
   # Change logging level in the script
   logging.basicConfig(level=logging.DEBUG)
   ```

2. **Disable Headless Mode:**
   ```python
   automator = SchwabLoginAutomator(headless=False)
   ```

3. **Add Breakpoints:**
   ```python
   import pdb; pdb.set_trace()  # Add where needed
   ```

## Customization

### Adding New Selectors
If Schwab updates their login form, you can add new selectors to the lists in `find_login_elements()`:

```python
username_selectors = [
    "input[name='LoginId']",
    "your_new_selector_here",
    # ... existing selectors
]
```

### Extending Functionality
You can extend the `SchwabLoginAutomator` class to add more features:

```python
def navigate_to_accounts(self):
    """Navigate to accounts page after login."""
    # Your implementation here
    pass
```

## Security Best Practices

1. **Never hardcode credentials** in the script
2. **Use environment variables** or secure credential managers
3. **Don't commit `.env` files** to version control
4. **Regularly update dependencies** for security patches
5. **Monitor script usage** and disable if suspicious activity occurs

## Legal Disclaimer

This script is provided for educational purposes only. Users are responsible for:
- Complying with Charles Schwab's terms of service
- Ensuring their use case is legal and ethical
- Protecting their account credentials
- Any consequences resulting from the use of this script

## Support

If you encounter issues:
1. Check the log files for error details
2. Review the screenshots for visual debugging
3. Ensure all dependencies are correctly installed
4. Verify your credentials and account status

## Dependencies

- `selenium>=4.15.0` - Web automation framework
- `webdriver-manager>=4.0.0` - Automatic WebDriver management
- `python-dotenv>=1.0.0` - Environment variable management