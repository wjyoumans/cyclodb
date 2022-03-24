Source code for building the [CycloDB](https://wjyoumans.github.io/cyclodb/) database and website.

To reproduce data for a specific field you can do (in `src` directory):

```sh
echo "cyclodata(63, 100)" | gp -D path=.:abelianbnf -D parisize=100000000 -q cyclodata.gp
```

This will output the raw data for the cyclotomic field of conductor 63 with a working precision of 100 decimal digits.
