#!/usr/bin/env python3
"""
Test script to verify Selenium and Firefox WebDriver setup
"""

from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from webdriver_manager.firefox import GeckoDriverManager
import time

def test_selenium_setup():
    """Test if Selenium can initialize Firefox WebDriver."""
    try:
        print("Testing Selenium and Firefox WebDriver setup...")
        
        # Firefox options for headless mode (since we're in a remote environment)
        firefox_options = Options()
        firefox_options.add_argument("--headless")
        firefox_options.add_argument("--no-sandbox")
        firefox_options.add_argument("--disable-dev-shm-usage")
        
        # Set up the service with WebDriver Manager
        print("Downloading and setting up GeckoDriver...")
        service = Service(GeckoDriverManager().install())
        
        # Initialize the driver
        print("Initializing Firefox WebDriver...")
        driver = webdriver.Firefox(service=service, options=firefox_options)
        
        # Test navigation
        print("Testing navigation to a test website...")
        driver.get("https://httpbin.org/html")
        
        # Get page title
        title = driver.title
        print(f"Page title: {title}")
        
        # Take a screenshot to verify it's working
        driver.save_screenshot("test_screenshot.png")
        print("Screenshot saved as test_screenshot.png")
        
        # Close the driver
        driver.quit()
        
        print("✅ Selenium setup test completed successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Error during setup test: {str(e)}")
        return False

if __name__ == "__main__":
    success = test_selenium_setup()
    if success:
        print("\n🎉 Your Schwab automation environment is ready!")
        print("You can now run: python schwab_login.py")
    else:
        print("\n⚠️  There were issues with the setup. Please check the error messages above.")