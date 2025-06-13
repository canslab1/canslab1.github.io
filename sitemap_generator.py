import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlparse
from datetime import datetime

BASE_URL = "https://canslab1.github.io/"
ALLOWED_EXTENSIONS = {'.html', '.pdf'}

def is_valid_file(url):
    parsed = urlparse(url)
    return any(parsed.path.endswith(ext) for ext in ALLOWED_EXTENSIONS)

def fetch_links(base_url):
    response = requests.get(base_url)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    links = set()

    for tag in soup.find_all("a", href=True):
        full_url = urljoin(base_url, tag["href"])
        if is_valid_file(full_url):
            links.add(full_url)
    return sorted(links)

def generate_sitemap_xml(urls, base_url):
    today = datetime.today().date()
    xml = ['<?xml version="1.0" encoding="UTF-8"?>']
    xml.append('<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">')
    
    for url in urls:
        xml.append("  <url>")
        xml.append(f"    <loc>{url}</loc>")
        xml.append(f"    <lastmod>{today}</lastmod>")
        xml.append("    <changefreq>daily</changefreq>")
        xml.append("    <priority>0.8</priority>")
        xml.append("  </url>")

    # Always include root
    xml.append("  <url>")
    xml.append(f"    <loc>{base_url}</loc>")
    xml.append(f"    <lastmod>{today}</lastmod>")
    xml.append("    <changefreq>daily</changefreq>")
    xml.append("    <priority>1.0</priority>")
    xml.append("  </url>")

    xml.append("</urlset>")
    return "\n".join(xml)

if __name__ == "__main__":
    try:
        urls = fetch_links(BASE_URL)
        sitemap_xml = generate_sitemap_xml(urls, BASE_URL)
        with open("sitemap.xml", "w", encoding="utf-8") as f:
            f.write(sitemap_xml)
        print("sitemap.xml generated successfully.")
    except Exception as e:
        print(f"Error occurred: {e}")
