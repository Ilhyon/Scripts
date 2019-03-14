# G4 Conserve

This project aims to predict G-quadruplexes (G4) in the whole human's transcriptome. We get our pG4, we analyze them.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To use all scripts in this project you will need to have python 2.7 and some libraries :
* pandas

You will also need to install G4RNA Screener from [Michelle Scott gitlab](http://gitlabscottgroup.med.usherbrooke.ca/J-Michel/g4rna_screener.git). Following is procedure to install it :

```
git clone http://gitlabscottgroup.med.usherbrooke.ca/J-Michel/g4rna_screener.git
sudo apt install python-pip
pip install biopython
pip install numpy
pip install pandas
pip install PyBrain
pip install regex
pip install scipy
cd PATH/TO/PYTHON/dist-packages/    or    cd PATH/TO/PYTHON/site-packages/
sudo -i
wget https://dev.mysql.com/get/Downloads/Connector-Python/mysql-connector-python-2.1.4.tar.gz
tar -xzf mysql-connector-python-2.1.4.tar.gz
cd mysql-connector-python-2.1.4
python setup.py install
```

### Installing

Except for the dependencies show earlier their is no installation. Their are only scripts.

Explain before G4 Screener
```
Give the example
```

And after

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc


