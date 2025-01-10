# DIYABC

We first of all need to install diyabc:

```bash
git clone --recurse-submodules https://github.com/diyabc/diyabc.git
cd diyabc
mkdir build
cd build
cmake ../
cmake --build . --config Release
```

You will find diyabc executable under `diyabc/build/src-JMC-C++/general`

Now we can just copy the binary file in our working directory:

```bash
cd ../..
cp diyabc/build/src-JMC-C++/general diyabc-RF-linux-v1.1.54
```

And we can optionally remove the diyabc repository we compiled as it it no longer necessary:

```bash
rm -rf diyabc/
```

And then we can test the installation by running:

```bash
./diyabc-RF-linux-v1.1.54 -h
```
