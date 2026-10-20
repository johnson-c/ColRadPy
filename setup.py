from distutils.core import setup
setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='colradpy',
    url='https://github.com/johnson-c/ColRadPy',
    author='Curt Johnson',
    author_email='caj0026@auburn.edu',
    # Needed to actually package something
    packages=['colradpy',],
    # Needed for dependencies
    install_requires=['numpy<2.5','scipy','matplotlib'],
    # *strongly* suggested for sharing
    version='1.1',
    description='python collisional radiative solver',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)
