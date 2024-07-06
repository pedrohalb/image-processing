# Código Base para processamento de imagem em C

## Pastas
1. **codigoBase** - contém um exemplo simples para calcular o negativo de uma imagem em formato pgm-ascii.
2. **utils** - contém a lib para ler e gravar imagens em formato pgm-ascii

## Instalação
1. Para o ambiente linux não é necessário instalar nenhum software adicional
2. Para windows:
   * Instalar o visualizador de imagens pgm-ascii [IrfanView](https://www.irfanview.com/64bit.htm);
   * Setar caminho para o visualizador
     - Sistema(Painel de Controle) -> Configurações Avançadas do Sistema -> Variáveis de Ambiente -> PATH (acrescentar o caminho "c:\Program Files\IrfanView")

## Compilar e Executar
1. Compilar a imagelib:
```
$ cd utils
$ make
```

2. Compilar e executar o exemplo
```
$ cd codigoBase
$ make
$ main img
```