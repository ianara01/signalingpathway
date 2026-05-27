# -*- coding: utf-8 -*-
"""
Created on Sun May 24 00:45:55 2026

@author: user
"""

import time
import requests

urls = [
    "https://reactome.org",
    "https://reactome.org/ContentService/",
    "https://reactome.org/ContentService/data/query/R-HSA-1234174",
]

for url in urls:
    print("\nURL:", url)
    try:
        r = requests.get(
            url,
            timeout=(10, 60),
            headers={"User-Agent": "S2SignalPathway/1.0"},
        )
        print("status:", r.status_code)
        print("text:", r.text[:120])
    except Exception as e:
        print(type(e).__name__, e)
    
    
"""
ReadTimeout HTTPSConnectionPool(host='reactome.org', port=443): Read timed out. (read timeout=10)
"""