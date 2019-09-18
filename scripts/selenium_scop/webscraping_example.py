from selenium import webdriver 
from selenium.webdriver.common.by import By 
from selenium.webdriver.support.ui import WebDriverWait 
from selenium.webdriver.support import expected_conditions as EC 
from selenium.common.exceptions import TimeoutException
import time


#option = webdriver.ChromeOptions()
#option.add_argument(" — incognito")

#browser = webdriver.Firefox(executable_path='/Users/admin/Applications/chromedriver', options=option)
browser = webdriver.Firefox(executable_path='/Users/admin/Applications/geckodriver')

browser.get('http://phosphatome.net/blast-scop/blast.html')

sequence = browser.find_element_by_name("SEQUENCE")
sequence.send_keys("MGAP")

sequence.submit()

time.sleep(5)

res = browser.find_element_by_xpath("(//pre)[3]")

print(res.text)

browser.quit() 



