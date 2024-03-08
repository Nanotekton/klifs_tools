import typer
from smiles_tools import modularize_ligand
from klifs_tools import check_allowed_connectivity
app = typer.Typer()

@app.command()
def test_ligand_blockization(smiles:str):
    '''Test "modularize_ligand" function'''
    result = modularize_ligand(smiles=smiles)
    print('blocks', '.'.join(result['blocks']))
    print('order', result['order'])
    for k,v in result['block_adj'].items():
        print(k, v)


@app.command()
def test_connectivity(test_sting:str):
    '''Test "check_allowed_connectivity" function'''
    keys = test_sting.split('-')
    print('Subpockets', KinFragLibAdj.keys())
    print('Connections:')#KinFragLibConn)
    for k in keys:
        print(k, KinFragLibAdj[k])
    print(check_allowed_connectivity(keys))
 

if __name__ == "__main__":
    app()
