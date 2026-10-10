"""Check local module syntax, document links, data paths, and generated-file inputs."""
from __future__ import annotations
import ast
import hashlib
import json
from pathlib import Path
import re
import sys
import unittest

ROOT=Path(__file__).resolve().parent
OPTIONAL={'flint','mpmath','numpy','scipy','matplotlib'}
RUNTIME={'.git','.venv','.lake','.deps','__pycache__','runs'}


def files():
    return [p for p in ROOT.rglob('*') if p.is_file() and p.relative_to(ROOT).parts[0] not in RUNTIME]


def walk_dicts(value):
    if isinstance(value,dict):
        yield value
        for x in value.values():yield from walk_dicts(x)
    elif isinstance(value,list):
        for x in value:yield from walk_dicts(x)


class RepositoryTests(unittest.TestCase):
    def test_all_python_syntax(self):
        for p in ROOT.glob('*.py'):
            with self.subTest(file=p.name):ast.parse(p.read_text(),filename=p.name)

    def test_local_imports_are_present(self):
        local={p.stem for p in ROOT.glob('*.py')}
        for p in ROOT.glob('*.py'):
            for node in ast.walk(ast.parse(p.read_text())):
                modules=[]
                if isinstance(node,ast.Import):modules=[a.name.split('.')[0] for a in node.names]
                if isinstance(node,ast.ImportFrom) and node.module:modules=[node.module.split('.')[0]]
                for name in modules:
                    with self.subTest(file=p.name,module=name):
                        self.assertIn(name,local|OPTIONAL|sys.stdlib_module_names|{'__future__'})

    def test_markdown_links_resolve(self):
        for p in ROOT.glob('*.md'):
            for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
                target=(p.parent/link.split('#')[0]).resolve()
                with self.subTest(file=p.name,link=link):
                    self.assertTrue(target.is_relative_to(ROOT))
                    self.assertTrue(target.is_file())

    def test_documented_modules_exist(self):
        for p in ROOT.glob('*.md'):
            for name in re.findall(r'\b([A-Za-z_]\w*\.(?:py|lean|md))\b',p.read_text()):
                with self.subTest(file=p.name,name=name):self.assertTrue((ROOT/name).is_file())

    def test_json_and_input_program_paths(self):
        program_count=0
        for p in files():
            if p.suffix!='.json':continue
            value=json.loads(p.read_text())
            for record in walk_dicts(value):
                if 'input_program' not in record:continue
                rel=record['input_program'];target=(ROOT/rel).resolve()
                with self.subTest(file=str(p.relative_to(ROOT)),program=rel):
                    self.assertTrue(target.is_relative_to(ROOT))
                    self.assertTrue(target.is_file())
                    self.assertEqual(target.suffix,'.spcf')
                program_count+=1
        self.assertEqual(program_count,78)

    def test_derivation_input_fingerprints(self):
        report=json.loads((ROOT/'results/derived/derivation_checks.json').read_text())
        for rel,expected in report['inputs_sha256'].items():
            with self.subTest(file=rel):
                self.assertEqual(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest(),expected)

    def test_formal_targets(self):
        targets=[]
        for p in ROOT.glob('Moment*.lean'):
            targets.extend(re.findall(r'^#print axioms (\S+)\s*$',p.read_text(),re.M))
        self.assertEqual(len(targets),18)
        self.assertEqual(len(set(targets)),18)

if __name__=='__main__':
    unittest.main(verbosity=2)
