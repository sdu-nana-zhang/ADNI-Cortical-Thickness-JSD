import sys,importlib.metadata
required=['numpy','pandas','openpyxl']
if '--figures' in sys.argv:required+=['scipy','matplotlib','nibabel','Pillow']
print(sys.version)
for name in required:
    try:print(name,importlib.metadata.version(name))
    except importlib.metadata.PackageNotFoundError:raise SystemExit('Missing Python package: '+name)
