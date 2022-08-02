import re

import selenium
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import os
import sapa as s

def cif_to_iso(cif_file,supercell=[2,0,0,0,2,0,0,0,2]):
    wd = os.getcwd()
    chrome_options = webdriver.ChromeOptions()
    prefs = {"download.default_directory": "%s" % wd}
    chrome_options.add_experimental_option("prefs", prefs)
    pathtofile = wd + "\\" + cif_file
    service = Service(executable_path=ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service,chrome_options=chrome_options)
    driver.get("https://stokes.byu.edu/iso/isodistort.php")
    driver.find_element(By.NAME, "toProcess").send_keys("%s"%pathtofile)
    driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
    driver.implicitly_wait(10)
    #Enter Supercell Coordinates
    supercell = [str(x) for x in supercell]
    basis = ["basis11","basis12","basis13","basis21","basis22","basis23","basis31","basis32","basis33"]
    for i in range(len(supercell)):
        driver.find_element(By.NAME, "%s"%basis[i]).clear()
        driver.find_element(By.NAME, "%s"%basis[i]).send_keys("%s"%supercell[i])

    #driver.find_element(By.NAME, "basis11").clear()
    #driver.find_element(By.NAME, "basis11").send_keys("2")
    #driver.find_element(By.NAME, "basis22").clear()
    #driver.find_element(By.NAME, "basis22").send_keys("2")
    #driver.find_element(By.NAME, "basis33").clear()
    #driver.find_element(By.NAME, "basis33").send_keys("2")
    #Pick symmetry
    select = Select(driver.find_element(By.NAME, "subgroupsym"))
    select.select_by_visible_text("1 P1 C1-1")
    driver.find_element(By.XPATH, "/html/body/div[2]/div[4]/form/h3/input").click()
    driver.implicitly_wait(10)
    window1 = driver.current_window_handle
    for window_handle in driver.window_handles:
        if window_handle != window1:
            driver.switch_to.window(window_handle)
            break
    driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
    driver.implicitly_wait(10)
    window2 = driver.current_window_handle
    for window_handle in driver.window_handles:
        if window_handle != window1 and window_handle != window2:
            driver.switch_to.window(window_handle)
            break
    #driver.find_element(By.XPATH,"/html/body/div[2]/form/input[46]").click()
    elements = driver.find_elements(By.CSS_SELECTOR, "input")
    for element in elements:
        if element.get_attribute("value") == "structurefile":
            element.click()
            break
    #driver.find_element(By.XPATH, "/html/body/div[2]/form/input[56]").click()
    driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
    time.sleep(5)
    #driver.find_element(By.XPATH, "/html/body/div[2]/form/input[44]").click()
    for element in elements:
        if element.get_attribute("value") == "isovizdistortion":
            element.click()
            break
    #driver.find_element(By.XPATH, "/html/body/div[2]/form/input[56]").click()
    driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
    time.sleep(5)
    driver.quit()
    fname = cif_file.split(".")
    isoname = fname[0] + "_iso.cif"
    os.rename("subgroup_cif.txt",isoname)

def test_opds(cif_file,kpoint,irrep,path):
        if not os.path.exists(f"./{path}"):
            os.makedirs(f"./{path}")
        wd = os.path.join(".",f"{path}")
        wd = os.path.abspath(wd)
        print(wd)
        chrome_options = webdriver.ChromeOptions()
        prefs = {"download.default_directory" : "{0}/".format(wd)}
        chrome_options.add_experimental_option("prefs", prefs)
        pathtofile = os.getcwd() + "\\" + cif_file
        service = Service(executable_path=ChromeDriverManager().install())          
        driver = webdriver.Chrome(service=service,chrome_options=chrome_options)
        driver.get("https://stokes.byu.edu/iso/isodistort.php")                     
        driver.find_element(By.NAME, "toProcess").send_keys("%s"%pathtofile)        
        driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()      
        driver.implicitly_wait(10)                                                  
        select = Select(driver.find_element(By.NAME, "kvec1"))
        #kpoint must be in vector form e.g. for X, kpoint = "(0,1/2,0)"
        for option in select.options:
            val = option.get_attribute("value")
            if kpoint in val:
                select.select_by_value(val)
                break
        window1 = driver.current_window_handle
        driver.find_element(By.XPATH,"/html/body/div[2]/div[3]/form[1]/h3/input").click()
        driver.implicitly_wait(10)
        for window_handle in driver.window_handles:
            if window_handle != window1:
                driver.switch_to.window(window_handle)
                break
        select_IR = Select(driver.find_element(By.NAME,"irrep1"))
        for option in select_IR.options:
            val = option.get_attribute("value")
            if irrep in val:
                select_IR.select_by_value(val)
                break
        #window2 = driver.current_window_handle
        #driver.find_element(By.XPATH,"/html/body/div[2]/form/input[39]").click()
        driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
        driver.implicitly_wait(10)
        #for window_handle in driver.window_handles:
        #    if window_handle != window1 and window_handle != window2:
        #        driver.switch_to.window(window_handle)
        #        break
        elements = driver.find_elements(By.CSS_SELECTOR,"input")
        opd_elements = []
        for element in elements:
            #print(element.get_attribute("type"))
            if element.get_attribute("type") == "radio":
                opd_elements.append(element)
        window3 = driver.current_window_handle
        for opd in opd_elements:
            opd.click()

            opdval = opd.get_attribute("value")
            regex = "\((.*?)\)"
            opdval = re.search(regex,opdval)
            match = opdval.group(1)
            print(match)
            #driver.find_element(By.XPATH,"/html/body/div[2]/form/input[46]").click()
            driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
            driver.implicitly_wait(10)  
            #print(opdval)
            for window_handle in driver.window_handles:
                if window_handle != window1 and window_handle != window3:
                    driver.switch_to.window(window_handle)
                    break
            #driver.find_element(By.XPATH,"/html/body/div[2]/form/input[51]")
            driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
            iso_options = driver.find_elements(By.CSS_SELECTOR,"input")
            for option in iso_options:
                if option.get_attribute("value") == "structurefile":
                    print(option.get_attribute("value"))
                    option.click()
                    
            #driver.find_element(By.XPATH,"/html/body/div[2]/form/input[61]").click()
            driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
            driver.implicitly_wait(10)

            time.sleep(5)
            os.rename(f"./{path}/subgroup_cif.txt",f"./{path}/{match}.cif")
            driver.close()
            driver.switch_to.window(window3)

        time.sleep(5)
        driver.implicitly_wait(10)

def findsym(cif_file,outfile):
    wd = os.getcwd()
    chrome_options = webdriver.ChromeOptions()
    prefs = {"download.default_directory": "%s" %wd, "download.prompt_for_download": False}
    chrome_options.add_experimental_option("prefs", prefs)
    pathtofile = wd + "\\" + cif_file
    service = Service(executable_path=ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, chrome_options=chrome_options)
    driver.get("https://stokes.byu.edu/iso/findsym.php")
    driver.find_element(By.NAME, "toProcess").send_keys("%s" % pathtofile)
    driver.find_element(By.XPATH, "//*[@class='btn btn-primary']").click()
    driver.implicitly_wait(10)
    driver.find_element(By.XPATH,"/html/body/form[1]/input[3]").click()
    time.sleep(5)
    os.rename("cif.txt",outfile)
