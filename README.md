
# PhyloDendron

Traditionally, phylogenetic analyses rely on command-line tools and scripting environments. [PhyloDendron.rocks](https://phylodendron.rocks) removes that barrier by providing an intuitive interface where users can upload data and visualize results directly in the browser.

A complete description can be found in the document Report/Final/Rapport_Baek_Slavova.pdf. 


 
## Authors

[@chelseabaek](https://www.github.com/chelseabaek) & [@mariyanas](https://www.github.com/mariyanas)

## Deployment

[PhyloDendron.rocks](https://phylodendron.rocks) is currently hosted on DigitalOcean.


## Run Locally

Clone the project

```bash
git clone https://github.com/chelseabaek/BIF7101-projet
```

Go to the project directory

```bash
cd app
```

Create a virtual environment

```bash
python -m venv venv
source venv/bin/activate
```

Install dependencies

```bash
pip install -r requirements.txt
```

Run with Docker

```bash
docker build -t phylodendron .
docker run -p 8080:8080 phylodendron
```

Then open

http://localhost:8080
