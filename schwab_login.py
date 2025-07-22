#!/usr/bin/env python3
"""
Schwab Login Automation Script using Selenium with Firefox
This script automates the login process for Charles Schwab's client portal.
"""

import os
import time
import logging
from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from webdriver_manager.firefox import GeckoDriverManager
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('schwab_login.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class SchwabLoginAutomator:
    """Automates login process for Charles Schwab client portal."""
    
    def __init__(self, headless=False, timeout=30):
        """
        Initialize the Schwab login automator.
        
        Args:
            headless (bool): Run browser in headless mode
            timeout (int): Default timeout for web elements
        """
        self.timeout = timeout
        self.driver = None
        self.wait = None
        self.headless = headless
        
    def setup_driver(self):
        """Set up Firefox WebDriver with optimized options."""
        try:
            # Firefox options
            firefox_options = Options()
            
            if self.headless:
                firefox_options.add_argument("--headless")
            
            # Security and performance options
            firefox_options.add_argument("--no-sandbox")
            firefox_options.add_argument("--disable-dev-shm-usage")
            firefox_options.add_argument("--disable-gpu")
            firefox_options.add_argument("--window-size=1920,1080")
            
            # Privacy settings
            firefox_options.set_preference("dom.webnotifications.enabled", False)
            firefox_options.set_preference("media.navigator.enabled", False)
            
            # Set up the service with WebDriver Manager
            service = Service(GeckoDriverManager().install())
            
            # Initialize the driver
            self.driver = webdriver.Firefox(service=service, options=firefox_options)
            self.driver.maximize_window()
            
            # Set up WebDriverWait
            self.wait = WebDriverWait(self.driver, self.timeout)
            
            logger.info("Firefox WebDriver initialized successfully")
            return True
            
        except Exception as e:
            logger.error(f"Failed to initialize WebDriver: {str(e)}")
            return False
    
    def navigate_to_login_page(self):
        """Navigate to Schwab login page."""
        try:
            login_url = "https://client.schwab.com/Areas/Access/Login"
            logger.info(f"Navigating to {login_url}")
            self.driver.get(login_url)
            
            # Wait for page to load
            self.wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
            
            # Additional wait for dynamic content
            time.sleep(3)
            
            logger.info("Successfully navigated to login page")
            return True
            
        except TimeoutException:
            logger.error("Timeout while loading login page")
            return False
        except Exception as e:
            logger.error(f"Error navigating to login page: {str(e)}")
            return False
    
    def find_login_elements(self):
        """Find and return login form elements."""
        try:
            # Common selectors for Schwab login form
            username_selectors = [
                "input[name='LoginId']",
                "input[id='LoginId']",
                "input[placeholder*='User']",
                "input[placeholder*='Login']",
                "input[type='text']",
                "#loginIdInput",
                ".login-input[name='LoginId']"
            ]
            
            password_selectors = [
                "input[name='Password']",
                "input[id='Password']",
                "input[type='password']",
                "#passwordInput",
                ".login-input[name='Password']"
            ]
            
            login_button_selectors = [
                "button[type='submit']",
                "input[type='submit']",
                "button[id*='login']",
                "button[class*='login']",
                ".login-button",
                "#loginButton"
            ]
            
            username_field = None
            password_field = None
            login_button = None
            
            # Try to find username field
            for selector in username_selectors:
                try:
                    username_field = self.driver.find_element(By.CSS_SELECTOR, selector)
                    logger.info(f"Found username field with selector: {selector}")
                    break
                except NoSuchElementException:
                    continue
            
            # Try to find password field
            for selector in password_selectors:
                try:
                    password_field = self.driver.find_element(By.CSS_SELECTOR, selector)
                    logger.info(f"Found password field with selector: {selector}")
                    break
                except NoSuchElementException:
                    continue
            
            # Try to find login button
            for selector in login_button_selectors:
                try:
                    login_button = self.driver.find_element(By.CSS_SELECTOR, selector)
                    logger.info(f"Found login button with selector: {selector}")
                    break
                except NoSuchElementException:
                    continue
            
            return username_field, password_field, login_button
            
        except Exception as e:
            logger.error(f"Error finding login elements: {str(e)}")
            return None, None, None
    
    def login(self, username, password):
        """
        Perform login with provided credentials.
        
        Args:
            username (str): Schwab username/login ID
            password (str): Schwab password
        """
        try:
            logger.info("Attempting to log in...")
            
            # Find login elements
            username_field, password_field, login_button = self.find_login_elements()
            
            if not all([username_field, password_field, login_button]):
                logger.error("Could not find all required login elements")
                return False
            
            # Clear and enter username
            username_field.clear()
            username_field.send_keys(username)
            logger.info("Username entered")
            
            # Clear and enter password
            password_field.clear()
            password_field.send_keys(password)
            logger.info("Password entered")
            
            # Small delay before clicking login
            time.sleep(1)
            
            # Click login button
            login_button.click()
            logger.info("Login button clicked")
            
            # Wait for navigation or error message
            time.sleep(5)
            
            # Check if login was successful
            current_url = self.driver.current_url
            if "login" not in current_url.lower() or "dashboard" in current_url.lower():
                logger.info("Login appears to be successful")
                return True
            else:
                logger.warning("Login may have failed - still on login page")
                return False
                
        except Exception as e:
            logger.error(f"Error during login: {str(e)}")
            return False
    
    def handle_two_factor_auth(self):
        """Handle two-factor authentication if required."""
        try:
            # Check for 2FA elements
            two_fa_selectors = [
                "input[name*='code']",
                "input[placeholder*='code']",
                "input[id*='verification']",
                ".verification-code",
                "#securityCode"
            ]
            
            for selector in two_fa_selectors:
                try:
                    code_field = self.driver.find_element(By.CSS_SELECTOR, selector)
                    logger.info("Two-factor authentication detected")
                    
                    # Prompt user for 2FA code
                    code = input("Enter your two-factor authentication code: ")
                    code_field.clear()
                    code_field.send_keys(code)
                    
                    # Find and click submit button
                    submit_button = self.driver.find_element(By.CSS_SELECTOR, "button[type='submit'], input[type='submit']")
                    submit_button.click()
                    
                    time.sleep(3)
                    return True
                    
                except NoSuchElementException:
                    continue
            
            return False
            
        except Exception as e:
            logger.error(f"Error handling 2FA: {str(e)}")
            return False
    
    def take_screenshot(self, filename="schwab_screenshot.png"):
        """Take a screenshot of the current page."""
        try:
            self.driver.save_screenshot(filename)
            logger.info(f"Screenshot saved as {filename}")
        except Exception as e:
            logger.error(f"Error taking screenshot: {str(e)}")
    
    def close(self):
        """Close the browser and clean up."""
        try:
            if self.driver:
                self.driver.quit()
                logger.info("Browser closed successfully")
        except Exception as e:
            logger.error(f"Error closing browser: {str(e)}")


def main():
    """Main function to run the Schwab login automation."""
    try:
        # Get credentials from environment variables or prompt user
        username = os.getenv('SCHWAB_USERNAME')
        password = os.getenv('SCHWAB_PASSWORD')
        
        if not username:
            username = input("Enter your Schwab username: ")
        
        if not password:
            import getpass
            password = getpass.getpass("Enter your Schwab password: ")
        
        # Initialize automator
        automator = SchwabLoginAutomator(headless=False)  # Set to True for headless mode
        
        # Set up driver
        if not automator.setup_driver():
            logger.error("Failed to set up WebDriver")
            return
        
        try:
            # Navigate to login page
            if not automator.navigate_to_login_page():
                logger.error("Failed to navigate to login page")
                return
            
            # Take screenshot of login page
            automator.take_screenshot("login_page.png")
            
            # Attempt login
            if automator.login(username, password):
                logger.info("Login successful!")
                
                # Handle 2FA if needed
                automator.handle_two_factor_auth()
                
                # Take screenshot after login
                automator.take_screenshot("after_login.png")
                
                # Keep browser open for a while to see the result
                logger.info("Keeping browser open for 30 seconds...")
                time.sleep(30)
                
            else:
                logger.error("Login failed")
                automator.take_screenshot("login_failed.png")
        
        finally:
            # Always close the browser
            automator.close()
    
    except KeyboardInterrupt:
        logger.info("Script interrupted by user")
    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")


if __name__ == "__main__":
    main()