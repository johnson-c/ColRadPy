from setuptools import setup
setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='colradpy',
    url='https://github.com/jladan/package_demo',
    author='Curt Johnson',
    author_email='caj0026@auburn.edu',
    # Needed to actually package something
    packages=['colradpy'],
    # Needed for dependencies
    install_requires=['numpy','collections','fractions','scipy','matplotlib','re','sys'],
    # *strongly* suggested for sharing
    version='1.1',
    description='python collisional radiative solver',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)
