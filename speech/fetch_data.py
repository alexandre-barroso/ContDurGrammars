"""Recover only the public annotation files required by this application."""
from pathlib import Path, PurePosixPath
import json
import tarfile
import urllib.request

HERE=Path(__file__).resolve().parent


def wanted(path):
    parts=PurePosixPath(path).parts
    return ((len(parts)==2 and parts[-1] in ['COPYING','README'])
            or (len(parts)==3 and parts[1]=='lab' and path.endswith('.lab'))
            or (len(parts)==4 and parts[1:3]==('festival','utts') and path.endswith('.utt'))
            or (len(parts)==3 and parts[1]=='etc' and parts[-1]=='txt.done.data'))


def main():
    sources=json.loads((HERE/'DATA_SOURCES.json').read_text())
    cache=HERE/'downloads';cache.mkdir(exist_ok=True)
    for source in sources['archives']:
        archive=cache/source['filename']
        if not archive.exists():
            urllib.request.urlretrieve(source['url'],archive)
        count=0
        with tarfile.open(archive,'r:bz2') as tar:
            for member in tar:
                parts=PurePosixPath(member.name).parts
                if '..' in parts or member.name.startswith('/'):raise ValueError('Unsafe archive path')
                if not member.isfile() or not wanted(member.name):continue
                destination=HERE/'data'/member.name
                destination.parent.mkdir(parents=True,exist_ok=True)
                destination.write_bytes(tar.extractfile(member).read());count+=1
        print(archive.name,count,'annotation and license files extracted')


if __name__=='__main__':main()
